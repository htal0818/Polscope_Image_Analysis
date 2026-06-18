function [BW, xc, yc, R_fit, polyXY, cache, geom] = segment_oocyte(Iseg, params, cache)
% SEGMENT_OOCYTE  Oocyte boundary segmentation (contour_retardance.m pipeline).
%
% Extracts the segmentation logic from contour_retardance.m into a reusable
% function for use by SCW_flows, boundary_flows, and other scripts.
%
% Inputs:
%   Iseg   - 2D image (double) for segmentation (e.g. 4-state sum)
%   params - struct with segmentation parameters:
%       .sigmaBlur            (default 1)    Gaussian blur sigma (px)
%       .closeRadius          (default 1)    morphological close disk radius (px)
%       .minArea              (default 5000) minimum object area (px^2)
%       .thresholdMode        (default 'adaptive') 'adaptive', 'edge', or 'gradient'
%       .segFromMask          (default true) true for 4-state (oocyte dark)
%       .adaptSensitivity     (default 0.7)  adaptthresh sensitivity [0..1]
%       .adaptNeighborhood    (default 201)  adaptthresh neighborhood (odd int, px)
%       .edgeMethod           (default 'Sobel') edge() method
%       .edgeDilateRadius     (default 2)    dilation radius (px) for edge gaps
%       .gradientPercentile   (default 70)   gradient magnitude percentile
%       .useCaching           (default true) enable mask caching between frames
%       .cacheIntensityThreshold (default 0.02) reuse mask if intensity change < this
%       .cacheForceRecalcEveryN  (default 25)   force recalc every N frames
%       .useActiveContour     (default false) optional Chan-Vese AC refinement
%       .activeContourIterations (default 80) AC iteration count
%       .activeContourMethod  (default 'Chan-Vese') activecontour method
%       .smoothMask           (default false) rebuild mask from smoothed boundary
%       .smoothBoundary       (default true) smooth returned boundary polygon
%       .nBoundaryPts         (default 720) number of arc-length boundary samples
%       .boundarySmoothFrac   (default 0.015) circular Gaussian window as contour fraction
%       .forceOuterEnvelope   (default false) rebuild a solid mask from radial outer envelope
%       .outerEnvelopePercentile (default 98) percentile radius per angular bin
%   cache  - struct for mask caching (pass [] on first call)
%
% Outputs:
%   BW     - binary mask of oocyte
%   xc, yc - circle fit center coordinates
%   R_fit  - circle fit radius (pixels)
%   polyXY - [N x 2] boundary polygon [x, y] from bwboundaries ([] if failed)
%   cache  - updated cache struct for next call
%   geom   - struct with smoothed-boundary geometry:
%            .polyXY, .arcLength_px, .curvature_pxInv, .radiusCurvature_px,
%            .theta_rad, .meanRadiusCurvature_px, .medianRadiusCurvature_px
%
% REQUIREMENTS:
%   - Image Processing Toolbox
%   - circfit.m

%% --- Default parameters ---
if nargin < 2 || isempty(params); params = struct(); end

sigmaBlur            = getfield_default(params, 'sigmaBlur', 10);
closeRadius          = getfield_default(params, 'closeRadius', 25);
minArea              = getfield_default(params, 'minArea', 5000);
thresholdMode        = getfield_default(params, 'thresholdMode', 'otsu');
segFromMask          = getfield_default(params, 'segFromMask', true);
adaptSensitivity     = getfield_default(params, 'adaptSensitivity', 0.7);
adaptNeighborhood    = getfield_default(params, 'adaptNeighborhood', 201);
edgeMethod           = getfield_default(params, 'edgeMethod', 'Canny');
edgeDilateRadius     = getfield_default(params, 'edgeDilateRadius', 1);
gradientPercentile   = getfield_default(params, 'gradientPercentile', 70);
useCaching           = getfield_default(params, 'useCaching', false);
cacheIntensityThreshold = getfield_default(params, 'cacheIntensityThreshold', 0.02);
cacheForceRecalcEveryN  = getfield_default(params, 'cacheForceRecalcEveryN', 25);
useActiveContour     = getfield_default(params, 'useActiveContour', false);
activeContourIterations = getfield_default(params, 'activeContourIterations', 80);
activeContourMethod  = getfield_default(params, 'activeContourMethod', 'edge');
smoothMask           = getfield_default(params, 'smoothMask', true);
smoothBoundary       = getfield_default(params, 'smoothBoundary', true);
nBoundaryPts         = getfield_default(params, 'nBoundaryPts', 720);
boundarySmoothFrac   = getfield_default(params, 'boundarySmoothFrac', 0.015);
forceOuterEnvelope   = getfield_default(params, 'forceOuterEnvelope', false);
outerEnvelopePercentile = getfield_default(params, 'outerEnvelopePercentile', 98);

%% --- Initialize cache ---
if isempty(cache)
    cache = struct();
    cache.prevBW      = [];
    cache.prevMeanInt = [];
    cache.frameNum    = 0;
    cache.hitCount    = 0;
    cache.recalcCount = 0;
end

cache.frameNum = cache.frameNum + 1;

%% --- Default outputs ---
BW = false(size(Iseg));
xc = NaN; yc = NaN; R_fit = NaN;
polyXY = [];
geom = empty_geom();

%% --- Decide whether to recalculate or reuse cached mask ---
needsRecalc = true;

if useCaching && ~isempty(cache.prevBW)
    if mod(cache.frameNum, cacheForceRecalcEveryN) == 1
        needsRecalc = true;   % forced drift correction
    else
        I_norm_check = Iseg / max(Iseg(:));
        meanIntCurrent = mean(I_norm_check(cache.prevBW), 'omitnan');
        intensityChange = abs(meanIntCurrent - cache.prevMeanInt) / (cache.prevMeanInt + eps);

        if intensityChange < cacheIntensityThreshold
            needsRecalc = false;  % cache hit
        end
    end
end

%% --- FULL MASK RECALCULATION or CACHE REUSE ---
if needsRecalc || ~useCaching
    I_blur = imgaussfilt(Iseg, sigmaBlur);
    I_norm = I_blur / max(I_blur(:));

    switch thresholdMode
        case 'adaptive'
            T = adaptthresh(Iseg);
            BW = imbinarize(Iseg, T);
            if segFromMask
                BW = ~BW;
            end

        case 'edge'
            edges = edge(Iseg, edgeMethod);
            se_edge = strel('disk', edgeDilateRadius);
            edges = imdilate(edges, se_edge);
            BW = imfill(edges, 'holes');
            if segFromMask && sum(BW(:)) > 0.5 * numel(BW)
                BW = ~BW;
            end
        case 'otsu'
            % Match contour_retardance.m: run Otsu on the blurred, normalized
            % image (I_norm), not raw Iseg. The blur washes out internal
            % oocyte texture so Otsu locks onto the gross egg shape, and
            % normalizing to [0,1] makes the threshold invariant to whether
            % Iseg is a 4-state sum or mean.
            Totsu = graythresh(I_norm);
            if segFromMask
                BW = I_norm < Totsu;
            else
                BW = I_norm > Totsu;
            end

        case 'gradient'
            [Gmag, ~] = imgradient(Iseg);
            thrG = prctile(Gmag(:), gradientPercentile);
            BW_edges = Gmag >= thrG;
            se_edge = strel('disk', edgeDilateRadius);
            BW_edges = imdilate(BW_edges, se_edge);
            BW = imfill(BW_edges, 'holes');
            if segFromMask && sum(BW(:)) > 0.5 * numel(BW)
                BW = ~BW;
            end

        otherwise
            error('segment_oocyte: Unknown thresholdMode: %s', thresholdMode);
    end

    % Aggressive morphological cleanup
    se = strel('disk', closeRadius);
    BW = imclose(BW, se);
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, minArea);

    % Fallback: gradient-based if threshold fails
    if ~any(BW(:))
        [Gmag, ~] = imgradient(Iseg);
        thrG = max(2*mean(Gmag(:)), prctile(Gmag(:), 90));
        BW = Gmag >= thrG;
        BW = imclose(BW, se);
        BW = imfill(BW, 'holes');
        BW = bwareaopen(BW, minArea);
    end

    % Keep largest connected component (fall back to cached mask on failure)
    L = bwlabel(BW, 8);
    if max(L(:)) >= 1
        S = regionprops(L, 'Area', 'Centroid');
        [~, iMax] = max([S.Area]);
        BW = (L == iMax);
    elseif ~isempty(cache.prevBW)
        BW = cache.prevBW;
    else
        % No boundary found
        return;
    end

    if useActiveContour
        BW_ac = activecontour(I_norm, BW, activeContourIterations, activeContourMethod);
        BW_ac = imclose(BW_ac, se);
        BW_ac = imfill(BW_ac, 'holes');
        BW_ac = bwareaopen(BW_ac, minArea);
        BW_ac = keep_largest_component(BW_ac);
        if any(BW_ac(:))
            BW = BW_ac;
        end
    end

    if forceOuterEnvelope
        BW_env = outer_envelope_mask(BW, nBoundaryPts, boundarySmoothFrac, outerEnvelopePercentile);
        BW_env = imfill(BW_env, 'holes');
        BW_env = bwareaopen(BW_env, minArea);
        BW_env = keep_largest_component(BW_env);
        if any(BW_env(:))
            BW = BW_env;
        end
    end

    if smoothMask || smoothBoundary
        [polySmooth, okSmooth] = smooth_mask_boundary(BW, nBoundaryPts, boundarySmoothFrac);
        if okSmooth && smoothMask
            [H, W] = size(BW);
            BW_smooth = poly2mask(polySmooth(:,1), polySmooth(:,2), H, W);
            BW_smooth = imfill(BW_smooth, 'holes');
            BW_smooth = bwareaopen(BW_smooth, minArea);
            BW_smooth = keep_largest_component(BW_smooth);
            if any(BW_smooth(:))
                BW = BW_smooth;
            end
        end
    end

    % Update cache
    if useCaching
        cache.prevBW = BW;
        I_norm_cache = Iseg / max(Iseg(:));
        cache.prevMeanInt = mean(I_norm_cache(BW), 'omitnan');
    end
    cache.recalcCount = cache.recalcCount + 1;

else
    % CACHE REUSE
    BW = cache.prevBW;
    cache.hitCount = cache.hitCount + 1;
end

%% --- Extract boundary contour ---
B = bwboundaries(BW);
if isempty(B)
    return;
end
[~, iLongest] = max(cellfun(@(p) size(p,1), B));
bnd = B{iLongest};
yb = bnd(:,1);
xb = bnd(:,2);

%% --- Circle fit for center & radius ---
[R_fit, xc, yc] = circfit(xb, yb);

%% --- Return smoothed boundary polygon and geometry if requested ---
if smoothBoundary
    [polyXY, okSmooth] = smooth_mask_boundary(BW, nBoundaryPts, boundarySmoothFrac);
    if ~okSmooth
        polyXY = [xb, yb];
    end
else
    polyXY = [xb, yb];
end

if size(polyXY,1) >= 5
    [R_fit, xc, yc] = circfit(polyXY(:,1), polyXY(:,2));
    geom = boundary_geometry(polyXY, xc, yc);
end

end

%% ========================== LOCAL HELPER ===================================
function v = getfield_default(s, fname, default)
    if isfield(s, fname)
        v = s.(fname);
    else
        v = default;
    end
end
function BWout = keep_largest_component(BWin)
    BWout = false(size(BWin));
    L = bwlabel(BWin, 8);
    if max(L(:)) < 1
        return;
    end
    S = regionprops(L, 'Area');
    [~, iMax] = max([S.Area]);
    BWout = (L == iMax);
end

function [polyXY, ok] = smooth_mask_boundary(BW, nPts, smoothFrac)
    polyXY = [];
    ok = false;

    B = bwboundaries(BW);
    if isempty(B)
        return;
    end

    [~, iLongest] = max(cellfun(@(p) size(p,1), B));
    bnd = B{iLongest};
    x = double(bnd(:,2));
    y = double(bnd(:,1));

    if numel(x) < 5
        return;
    end

    if x(1) ~= x(end) || y(1) ~= y(end)
        x(end+1) = x(1);
        y(end+1) = y(1);
    end

    ds = hypot(diff(x), diff(y));
    keep = [true; ds > 0];
    x = x(keep);
    y = y(keep);
    ds = hypot(diff(x), diff(y));
    cumLen = [0; cumsum(ds)];
    totalLen = cumLen(end);

    if totalLen <= 0 || numel(cumLen) < 5
        return;
    end

    nPts = max(64, round(nPts));
    sUniform = linspace(0, totalLen, nPts + 1)';
    sUniform(end) = [];

    xUniform = interp1(cumLen, x, sUniform, 'pchip');
    yUniform = interp1(cumLen, y, sUniform, 'pchip');

    smoothW = max(5, round(nPts * smoothFrac));
    if mod(smoothW, 2) == 0
        smoothW = smoothW + 1;
    end
    smoothW = min(smoothW, floor((nPts - 1) / 2));

    if smoothW >= 5
        xPad = [xUniform(end-smoothW+1:end); xUniform; xUniform(1:smoothW)];
        yPad = [yUniform(end-smoothW+1:end); yUniform; yUniform(1:smoothW)];
        xPad = smoothdata(xPad, 'gaussian', smoothW);
        yPad = smoothdata(yPad, 'gaussian', smoothW);
        xUniform = xPad(smoothW+1:smoothW+nPts);
        yUniform = yPad(smoothW+1:smoothW+nPts);
    end

    [H, W] = size(BW);
    xUniform = min(max(xUniform, 1), W);
    yUniform = min(max(yUniform, 1), H);
    polyXY = [xUniform(:), yUniform(:)];
    ok = true;
end

function BWenv = outer_envelope_mask(BW, nPts, smoothFrac, radiusPercentile)
    BWenv = false(size(BW));
    BW = keep_largest_component(BW);
    if ~any(BW(:))
        return;
    end

    S = regionprops(BW, 'Centroid');
    if isempty(S)
        return;
    end
    c = S(1).Centroid;
    xc = c(1);
    yc = c(2);

    [yy, xx] = find(BW);
    x = double(xx);
    y = double(yy);
    theta = atan2(y - yc, x - xc);
    theta(theta < 0) = theta(theta < 0) + 2*pi;
    r = hypot(x - xc, y - yc);

    nPts = max(64, round(nPts));
    edges = linspace(0, 2*pi, nPts + 1);
    rEnv = nan(nPts, 1);
    for ii = 1:nPts
        inBin = theta >= edges(ii) & theta < edges(ii+1);
        if any(inBin)
            rEnv(ii) = prctile(r(inBin), radiusPercentile);
        end
    end

    good = isfinite(rEnv);
    if nnz(good) < 8
        return;
    end

    idx = (1:nPts)';
    idxGood = idx(good);
    rGood = rEnv(good);
    idxExt = [idxGood - nPts; idxGood; idxGood + nPts];
    rExt = [rGood; rGood; rGood];
    rEnv = interp1(idxExt, rExt, idx, 'pchip');

    smoothW = max(5, round(nPts * smoothFrac));
    if mod(smoothW, 2) == 0
        smoothW = smoothW + 1;
    end
    smoothW = min(smoothW, floor((nPts - 1) / 2));
    if smoothW >= 5
        rPad = [rEnv(end-smoothW+1:end); rEnv; rEnv(1:smoothW)];
        rPad = smoothdata(rPad, 'gaussian', smoothW);
        rEnv = rPad(smoothW+1:smoothW+nPts);
    end

    th = linspace(0, 2*pi, nPts + 1)';
    th(end) = [];
    xPoly = xc + rEnv .* cos(th);
    yPoly = yc + rEnv .* sin(th);

    [H, W] = size(BW);
    xPoly = min(max(xPoly, 1), W);
    yPoly = min(max(yPoly, 1), H);
    BWenv = poly2mask(xPoly, yPoly, H, W);
end

function geom = boundary_geometry(polyXY, xc, yc)
    geom = empty_geom();
    n = size(polyXY, 1);
    if n < 5
        return;
    end

    x = polyXY(:,1);
    y = polyXY(:,2);
    xNext = circshift(x, -1);
    yNext = circshift(y, -1);
    dsEach = hypot(xNext - x, yNext - y);
    ds_px = median(dsEach(dsEach > 0), 'omitnan');
    if ~isfinite(ds_px) || ds_px <= 0
        return;
    end

    xPrev = circshift(x, 1);
    yPrev = circshift(y, 1);
    dx_ds = (xNext - xPrev) / (2 * ds_px);
    dy_ds = (yNext - yPrev) / (2 * ds_px);
    d2x_ds2 = (xNext - 2*x + xPrev) / (ds_px^2);
    d2y_ds2 = (yNext - 2*y + yPrev) / (ds_px^2);

    denom = (dx_ds.^2 + dy_ds.^2).^(3/2);
    kappa = (dx_ds .* d2y_ds2 - dy_ds .* d2x_ds2) ./ (denom + eps);
    radiusCurv = 1 ./ (abs(kappa) + eps);

    theta = atan2(y - yc, x - xc);
    theta(theta < 0) = theta(theta < 0) + 2*pi;
    arcLength = (0:n-1)' * ds_px;

    geom.polyXY = polyXY;
    geom.arcLength_px = arcLength;
    geom.curvature_pxInv = kappa;
    geom.radiusCurvature_px = radiusCurv;
    geom.theta_rad = theta;
    geom.meanRadiusCurvature_px = mean(radiusCurv, 'omitnan');
    geom.medianRadiusCurvature_px = median(radiusCurv, 'omitnan');
end

function geom = empty_geom()
    geom = struct( ...
        'polyXY', [], ...
        'arcLength_px', [], ...
        'curvature_pxInv', [], ...
        'radiusCurvature_px', [], ...
        'theta_rad', [], ...
        'meanRadiusCurvature_px', NaN, ...
        'medianRadiusCurvature_px', NaN);
end

function result = measure_contour_retardance_test(Iraw, opts)
% MEASURE_CONTOUR_RETARDANCE_TEST  Test version with peak-along-inward-normal cortex measurement.
%
%   result = measure_contour_retardance_test(Iraw, opts)
%
%   Same segmentation pipeline as the original measure_contour_retardance, but
%   the cortex retardance at each angular position is the PEAK retardance along
%   the inward normal at that boundary point (within peakSearchDepth_um), not a
%   sample at a fixed inward offset. The "profile with the highest retardance
%   values is representative of the cortex" — this function takes the per-angle
%   peak of every inward radial profile.
%
%   This is a sandbox / test version intended for use by batch_distribution_test.m.
%   The working contour_retardance.m and any production batch script are not
%   touched.
%
%   INPUTS
%     Iraw  — [H x W] double, raw pixel values (e.g. 0..65535 for 16-bit)
%     opts  — struct with fields (all optional, defaults shown):
%               retardance_ceiling_nm  (50)     Polscope ceiling in nm
%               bit_depth              (16)     image bit depth
%               sigmaBlur              (20)     Gaussian blur sigma for segmentation
%               closeRadius            (25)     morphological close disk radius
%               minArea                (5000)   minimum object area in px^2
%               peakSearchDepth_um     (5)      max depth (um) along inward normal
%                                               to search for cortical peak
%               nBoundaryPts           (500)    uniformly spaced boundary points
%               maxDepth_um            (15)     how far inward to scan (um)
%               depthStep_um           (0.2)    step size along normal (um)
%               px_per_um              (6.25)   spatial calibration
%               thresholdMode          ('otsu') 'otsu', 'fixed', 'percentile',
%                                               or 'adaptive'
%               fixedThreshold         (500)    for 'fixed' mode
%               percentileThreshold    (30)     for 'percentile' mode
%               adaptiveSensitivity    (0.5)    for 'adaptive' mode (0..1)
%               prevBW                 ([])     previous frame mask for fallback
%               Iseg                   ([])     separate image for segmentation
%                                               (e.g. avg of State1-4). If empty,
%                                               Iraw is used for segmentation.
%
%   OUTPUT
%     result — struct with fields:
%               contourValues           — [nBoundaryPts x 1] per-angle PEAK
%                                          retardance (nm) along inward normal
%               contourMean/Std/Max/Min — scalars over contourValues
%               peakDepth_um            — [nBoundaryPts x 1] depth of each peak
%               normalProfilesPerAngle  — [nBoundaryPts x nDepth] full profiles
%               depthAxis_um            — [1 x nDepth]
%               polyX, polyY            — smooth boundary points
%               nx, ny                  — inward normals
%               xb_orig, yb_orig        — raw boundary contour (unsmoothed)
%               xc, yc                  — circle-fit center
%               R_fit                   — circle-fit radius (px)
%               BW                      — binary mask of the oocyte
%               Iret                    — retardance image in nm
%               success                 — logical
%
%   Requires: circfit.m (in this repository), Image Processing Toolbox

    %% Defaults
    if nargin < 2; opts = struct(); end
    def = struct( ...
        'retardance_ceiling_nm', 50, ...
        'bit_depth',             16, ...
        'sigmaBlur',             20, ...
        'closeRadius',           25, ...
        'minArea',               5000, ...
        'peakSearchDepth_um',    5, ...
        'nBoundaryPts',          500, ...
        'maxDepth_um',           15, ...
        'depthStep_um',          0.2, ...
        'px_per_um',             6.25, ...
        'thresholdMode',         'otsu', ...
        'fixedThreshold',        500, ...
        'percentileThreshold',   30, ...
        'adaptiveSensitivity',   0.5, ...
        'prevBW',                [], ...
        'Iseg',                  []);
    flds = fieldnames(def);
    for k = 1:numel(flds)
        if ~isfield(opts, flds{k})
            opts.(flds{k}) = def.(flds{k});
        end
    end

    %% Convert to retardance (nm)
    maxPixVal = 2^opts.bit_depth - 1;
    Iret = (Iraw / maxPixVal) * opts.retardance_ceiling_nm;
    [H, W] = size(Iraw);

    %% Depth axis along inward normals
    um_per_px    = 1 / opts.px_per_um;
    depthAxis_um = 0 : opts.depthStep_um : opts.maxDepth_um;
    nDepth       = numel(depthAxis_um);
    depthAxis_px = depthAxis_um / um_per_px;

    %% Boundary detection (segmentation source: Iseg if provided, else Iraw)
    if ~isempty(opts.Iseg)
        Iseg = opts.Iseg;
    else
        Iseg = Iraw;
    end

    switch opts.thresholdMode
        case 'otsu'
            I_blur = imgaussfilt(Iseg, opts.sigmaBlur);
            I_norm = I_blur / max(I_blur(:));
            BW = I_norm > graythresh(I_norm);
        case 'fixed'
            I_blur = imgaussfilt(Iseg, opts.sigmaBlur);
            BW = I_blur > opts.fixedThreshold;
        case 'percentile'
            I_blur = imgaussfilt(Iseg, opts.sigmaBlur);
            pVal = prctile(I_blur(:), opts.percentileThreshold);
            BW = I_blur > pVal;
        case 'adaptive'
            I_blur = imgaussfilt(Iseg, opts.sigmaBlur);
            I_norm = I_blur / max(I_blur(:));
            T = adaptthresh(I_norm, opts.adaptiveSensitivity);
            BW = imbinarize(I_norm, T);
        otherwise
            error('Unknown thresholdMode: %s', opts.thresholdMode);
    end

    % Morphological cleanup
    se = strel('disk', opts.closeRadius);
    BW = imclose(BW, se);
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, opts.minArea);

    % Fallback: gradient-based if threshold collapses to empty
    if ~any(BW(:))
        [Gmag, ~] = imgradient(I_blur);
        thrG = max(2*mean(Gmag(:)), prctile(Gmag(:), 80));
        BW = Gmag >= thrG;
        BW = imclose(BW, se);
        BW = imfill(BW, 'holes');
        BW = bwareaopen(BW, opts.minArea);
    end

    % Keep largest connected component
    L = bwlabel(BW, 8);
    if max(L(:)) >= 1
        S = regionprops(L, 'Area');
        [~, iMax] = max([S.Area]);
        BW = (L == iMax);
    elseif ~isempty(opts.prevBW)
        BW = opts.prevBW;
    else
        result = make_empty_result(Iret, BW, depthAxis_um);
        return;
    end

    %% Extract boundary contour
    B = bwboundaries(BW);
    if isempty(B)
        result = make_empty_result(Iret, BW, depthAxis_um);
        return;
    end
    [~, iLongest] = max(cellfun(@(p) size(p,1), B));
    bnd     = B{iLongest};
    yb_orig = bnd(:,1);
    xb_orig = bnd(:,2);

    %% Circle fit for center & radius
    [R_fit, xc, yc] = circfit(xb_orig, yb_orig);

    %% Smooth boundary contour (uniform arc-length resample + circular Gaussian smooth)
    cumLen   = [0; cumsum(sqrt(diff(xb_orig).^2 + diff(yb_orig).^2))];
    uniformS = linspace(0, cumLen(end), opts.nBoundaryPts + 1)';
    uniformS(end) = [];
    polyX = interp1(cumLen, xb_orig, uniformS, 'pchip')';
    polyY = interp1(cumLen, yb_orig, uniformS, 'pchip')';

    smoothW = max(5, round(opts.nBoundaryPts * 0.10));
    xPad    = [polyX(end-smoothW+1:end), polyX, polyX(1:smoothW)];
    yPad    = [polyY(end-smoothW+1:end), polyY, polyY(1:smoothW)];
    xPad    = smoothdata(xPad, 'gaussian', smoothW);
    yPad    = smoothdata(yPad, 'gaussian', smoothW);
    polyX   = xPad(smoothW+1 : smoothW+opts.nBoundaryPts);
    polyY   = yPad(smoothW+1 : smoothW+opts.nBoundaryPts);

    %% Analytic normals from distance transform gradient
    per      = bwperim(BW);
    D_full   = bwdist(per);
    [Gy, Gx] = imgradientxy(D_full, 'central');

    F_Gx = griddedInterpolant({1:H, 1:W}, Gx, 'linear', 'nearest');
    F_Gy = griddedInterpolant({1:H, 1:W}, Gy, 'linear', 'nearest');
    nx = -F_Gx(polyY, polyX);   % inward = negative gradient
    ny = -F_Gy(polyY, polyX);
    nmag = sqrt(nx.^2 + ny.^2) + eps;
    nx = nx ./ nmag;
    ny = ny ./ nmag;

    % Verify inward orientation (dot with center direction)
    toCenter_x = xc - polyX;  toCenter_y = yc - polyY;
    dot_check  = nx .* toCenter_x + ny .* toCenter_y;
    nx(dot_check < 0) = -nx(dot_check < 0);
    ny(dot_check < 0) = -ny(dot_check < 0);

    %% Per-angle peak retardance along inward normals
    F = griddedInterpolant({1:H, 1:W}, Iret, 'linear', 'nearest');

    peakSearchIdx = find(depthAxis_um <= opts.peakSearchDepth_um);
    contourValues = nan(opts.nBoundaryPts, 1);
    peakDepth_um  = nan(opts.nBoundaryPts, 1);
    normalProfilesPerAngle = nan(opts.nBoundaryPts, nDepth);

    for bi = 1:opts.nBoundaryPts
        xs_n = polyX(bi) + depthAxis_px * nx(bi);
        ys_n = polyY(bi) + depthAxis_px * ny(bi);

        inBounds = xs_n >= 1 & xs_n <= W & ys_n >= 1 & ys_n <= H;
        if ~any(inBounds); continue; end

        vals = F(ys_n(inBounds), xs_n(inBounds));
        normalProfilesPerAngle(bi, inBounds) = vals(:)';

        searchIdx = intersect(peakSearchIdx, find(inBounds));
        if ~isempty(searchIdx)
            searchVals = F(ys_n(searchIdx), xs_n(searchIdx));
            [contourValues(bi), pidx] = max(searchVals);
            peakDepth_um(bi) = depthAxis_um(searchIdx(pidx));
        end
    end

    %% Pack result
    result.contourValues          = contourValues;
    result.contourMean            = mean(contourValues, 'omitnan');
    result.contourStd             = std(contourValues, 0, 'omitnan');
    result.contourMax             = max(contourValues);
    result.contourMin             = min(contourValues);
    result.peakDepth_um           = peakDepth_um;
    result.normalProfilesPerAngle = normalProfilesPerAngle;
    result.depthAxis_um           = depthAxis_um;
    result.polyX                  = polyX;
    result.polyY                  = polyY;
    result.nx                     = nx;
    result.ny                     = ny;
    result.xb_orig                = xb_orig;
    result.yb_orig                = yb_orig;
    result.xc                     = xc;
    result.yc                     = yc;
    result.R_fit                  = R_fit;
    result.BW                     = BW;
    result.Iret                   = Iret;
    result.success                = true;
end

function result = make_empty_result(Iret, BW, depthAxis_um)
    result.contourValues          = [];
    result.contourMean            = NaN;
    result.contourStd             = NaN;
    result.contourMax             = NaN;
    result.contourMin             = NaN;
    result.peakDepth_um           = [];
    result.normalProfilesPerAngle = [];
    result.depthAxis_um           = depthAxis_um;
    result.polyX                  = [];
    result.polyY                  = [];
    result.nx                     = [];
    result.ny                     = [];
    result.xb_orig                = [];
    result.yb_orig                = [];
    result.xc                     = NaN;
    result.yc                     = NaN;
    result.R_fit                  = NaN;
    result.BW                     = BW;
    result.Iret                   = Iret;
    result.success                = false;
end

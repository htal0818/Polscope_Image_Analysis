function [BW, xc, yc, R_fit, polyXY, cache] = segment_oocyte(Iseg, params, cache)
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
%   cache  - struct for mask caching (pass [] on first call)
%
% Outputs:
%   BW     - binary mask of oocyte
%   xc, yc - circle fit center coordinates
%   R_fit  - circle fit radius (pixels)
%   polyXY - [N x 2] boundary polygon [x, y] from bwboundaries ([] if failed)
%   cache  - updated cache struct for next call
%
% REQUIREMENTS:
%   - Image Processing Toolbox
%   - circfit.m

%% --- Default parameters ---
if nargin < 2 || isempty(params); params = struct(); end

sigmaBlur            = getfield_default(params, 'sigmaBlur', 1);
closeRadius          = getfield_default(params, 'closeRadius', 1);
minArea              = getfield_default(params, 'minArea', 5000);
thresholdMode        = getfield_default(params, 'thresholdMode', 'adaptive');
segFromMask          = getfield_default(params, 'segFromMask', true);
adaptSensitivity     = getfield_default(params, 'adaptSensitivity', 0.7);
adaptNeighborhood    = getfield_default(params, 'adaptNeighborhood', 201);
edgeMethod           = getfield_default(params, 'edgeMethod', 'Sobel');
edgeDilateRadius     = getfield_default(params, 'edgeDilateRadius', 2);
gradientPercentile   = getfield_default(params, 'gradientPercentile', 70);
useCaching           = getfield_default(params, 'useCaching', true);
cacheIntensityThreshold = getfield_default(params, 'cacheIntensityThreshold', 0.02);
cacheForceRecalcEveryN  = getfield_default(params, 'cacheForceRecalcEveryN', 25);

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
            edges = imerode(edges, se_edge);
            BW = imfill(edges, 'holes');
            if segFromMask && sum(BW(:)) > 0.5 * numel(BW)
                BW = ~BW;
            end

        case 'gradient'
            [Gmag, ~] = imgradient(Iseg);
            thrG = prctile(Gmag(:), gradientPercentile);
            BW_edges = Gmag >= thrG;
            se_edge = strel('disk', edgeDilateRadius);
            BW_edges = imerode(BW_edges, se_edge);
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

    % Reject components touching the image border (vignette / dark frame
    % corners / acquisition-mask artifacts). The oocyte is centrally
    % located, so border-touching blobs are never the object of interest.
    % Guarded: if clearing the border kills everything, keep the original.
    BW_cleared = imclearborder(BW);
    if any(BW_cleared(:))
        BW = BW_cleared;
    end

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

%% --- Return boundary polygon as [x, y] ---
polyXY = [xb, yb];

end

%% ========================== LOCAL HELPER ===================================
function v = getfield_default(s, fname, default)
    if isfield(s, fname)
        v = s.(fname);
    else
        v = default;
    end
end

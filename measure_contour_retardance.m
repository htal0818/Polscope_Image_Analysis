function result = measure_contour_retardance(Iraw, opts)
% MEASURE_CONTOUR_RETARDANCE  Segment an oocyte and measure retardance at its contour.
%
%   result = measure_contour_retardance(Iraw, opts)
%
%   Segmentation uses Laplacian-of-Gaussian ridge detection to build an
%   initial seed, refined by an edge-based active contour.
%
%   INPUTS
%     Iraw  — [H x W] double, raw pixel values (e.g. 0..65535 for 16-bit)
%     opts  — struct with fields (all optional, defaults shown):
%               retardance_ceiling_nm    (50)       Polscope ceiling in nm
%               bit_depth                (16)       image bit depth
%               smoothSigma              (1.8)      Gaussian blur sigma for segmentation
%               ridgeSigma               (2.0)      LoG kernel sigma for ridge detection
%               closeRadius              (12)       morphological close disk radius for seed
%               seedDilateRadius         (3)        dilate seed before active contour
%               minArea                  (5000)     minimum object area in px^2
%               boundaryInset_px         (1)        shift boundary inward onto cortex
%               activeContourIterations  (200)      active contour iterations
%               activeContourSmoothness  (1.5)      active contour smooth factor
%               edgeContractionBias      (0.0)      active contour contraction bias
%               um_per_px                (1)        spatial calibration
%               profileMaxDepth_um       (50)       max inward depth for radial profile
%               profileDepthStep_um      (1)        inward radial bin size
%               peakSearchDepthRange_um  ([0 50])   search range for peak ring
%               prevBW                   ([])       previous frame mask for fallback
%               Iseg                     ([])       separate image for segmentation;
%                                                   if empty, Iraw is used
%
%   OUTPUT
%     result — struct with fields:
%               contourValues, contourMean, contourStd, contourMax, contourMin
%               xb, yb          — boundary point coordinates (after inset)
%               xc, yc          — circle-fit center
%               R_fit           — circle-fit radius (px)
%               BW              — final binary mask
%               BW_initial      — pre-active-contour seed mask
%               ridgeResponse   — LoG ridge response image
%               gradMask        — thresholded ridge mask
%               edgeMask        — Canny fallback edge mask
%               boundarySupport — (unused, for compatibility)
%               peakRingValues, peakRingDepth_um, peakRingMean, etc.
%               Iret            — retardance image in nm
%               success         — logical
%
%   Requires: circfit.m, Image Processing Toolbox

    %% Defaults
    if nargin < 2; opts = struct(); end
    def = struct( ...
        'retardance_ceiling_nm', 50, ...
        'bit_depth',             16, ...
        'smoothSigma',           1.8, ...
        'ridgeSigma',            2.0, ...
        'closeRadius',           12, ...
        'seedDilateRadius',      3, ...
        'minArea',               5000, ...
        'boundaryInset_px',      1, ...
        'activeContourIterations', 200, ...
        'activeContourSmoothness', 1.5, ...
        'edgeContractionBias',   0.0, ...
        'um_per_px',             1, ...
        'profileMaxDepth_um',    50, ...
        'profileDepthStep_um',   1, ...
        'peakSearchDepthRange_um', [0 50], ...
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

    %% Select segmentation image
    if ~isempty(opts.Iseg)
        Iseg = opts.Iseg;
    else
        Iseg = Iraw;
    end

    %% LoG ridge-based segmentation
    finitePixels = Iseg(isfinite(Iseg));
    displayLimits = prctile(finitePixels, [0.1 99.9]);
    Inorm = (Iseg - displayLimits(1)) / (displayLimits(2) - displayLimits(1) + eps);
    Inorm = min(max(Inorm, 0), 1);
    segImg = imgaussfilt(Inorm, opts.smoothSigma);

    kernelSize = 2*ceil(3*opts.ridgeSigma) + 1;
    logKernel = fspecial('log', kernelSize, opts.ridgeSigma);
    signedLoG = imfilter(segImg, logKernel, 'replicate', 'conv');
    ridgeResponse = max(-signedLoG, 0);
    ridgeScale = prctile(ridgeResponse(:), 99.8);
    ridgeResponse = min(ridgeResponse / (ridgeScale + eps), 1);

    borderWidth = max(8, round(0.05 * min(size(Iseg))));
    borderMask = false(size(Iseg));
    borderMask([1:borderWidth, end-borderWidth+1:end], :) = true;
    borderMask(:, [1:borderWidth, end-borderWidth+1:end]) = true;
    background = ridgeResponse(borderMask);
    noiseThreshold = median(background) + 6 * 1.4826 * mad(background, 1);
    ridgeThreshold = min(max(graythresh(ridgeResponse), noiseThreshold), 0.98);

    ridgeMask = ridgeResponse > ridgeThreshold;
    gradMask = ridgeMask;
    edgeMask = false(size(Iseg));
    boundarySupport = false(size(Iseg));

    seed = make_filled_seed(ridgeMask, opts.closeRadius);

    % Canny fallback
    if isempty(seed)
        edgeMask = edge(segImg, 'Canny', [], max(1, opts.smoothSigma));
        seed = make_filled_seed(edgeMask, opts.closeRadius);
    end

    if isempty(seed)
        if ~isempty(opts.prevBW)
            BW = opts.prevBW;
            BW_initial = BW;
        else
            result = make_empty_result(Iret, false(H, W), gradMask, edgeMask, boundarySupport);
            result.ridgeResponse = ridgeResponse;
            return;
        end
    else
        if opts.seedDilateRadius > 0
            seed = imdilate(seed, strel('disk', opts.seedDilateRadius, 0));
        end
        seed = imfill(seed, 'holes');
        BW_initial = seed;

        BW = activecontour(segImg, seed, opts.activeContourIterations, 'edge', ...
            'SmoothFactor', opts.activeContourSmoothness, ...
            'ContractionBias', opts.edgeContractionBias);
        BW = imfill(BW, 'holes');
        BW = bwareaopen(BW, opts.minArea);

        if any(BW(:))
            BW = keep_largest_component(BW);
        elseif ~isempty(opts.prevBW)
            BW = opts.prevBW;
        else
            result = make_empty_result(Iret, BW_initial, gradMask, edgeMask, boundarySupport);
            result.ridgeResponse = ridgeResponse;
            return;
        end
    end

    %% Extract boundary contour
    B = bwboundaries(BW);
    if isempty(B)
        result = make_empty_result(Iret, BW, gradMask, edgeMask, boundarySupport);
        result.ridgeResponse = ridgeResponse;
        return;
    end
    [~, iLongest] = max(cellfun(@(p) size(p,1), B));
    bnd = B{iLongest};
    yb = bnd(:,1);
    xb = bnd(:,2);

    %% Circle fit for center & radius
    [R_fit, xc, yc] = circfit(xb, yb);

    % Shrink boundary inward onto cortical ring center
    dx = xb - xc;  dy = yb - yc;
    dist = sqrt(dx.^2 + dy.^2);
    shrink = max(dist - opts.boundaryInset_px, 1) ./ dist;
    xb = xc + dx .* shrink;
    yb = yc + dy .* shrink;

    %% Interpolate retardance (nm) at boundary pixel locations
    F = griddedInterpolant({1:H, 1:W}, Iret, 'linear', 'nearest');
    ib = F(yb, xb);

    %% Peak radial ring statistics
    peakRing = compute_peak_ring_stats(Iret, BW, opts);

    %% Pack result
    result.contourValues = ib;
    result.contourMean   = mean(ib, 'omitnan');
    result.contourStd    = std(ib, 'omitnan');
    result.contourMax    = max(ib);
    result.contourMin    = min(ib);
    result.xb            = xb;
    result.yb            = yb;
    result.xc            = xc;
    result.yc            = yc;
    result.R_fit         = R_fit;
    result.BW            = BW;
    result.BW_initial    = BW_initial;
    result.ridgeResponse = ridgeResponse;
    result.gradMask      = gradMask;
    result.edgeMask      = edgeMask;
    result.boundarySupport = boundarySupport;
    result.depthAxis_um    = peakRing.depthAxis_um;
    result.distProfile_nm  = peakRing.distProfile_nm;
    result.peakRingDepth_um = peakRing.depth_um;
    result.peakRingValues  = peakRing.values_nm;
    result.peakRingMean    = peakRing.mean_nm;
    result.peakRingMedian  = peakRing.median_nm;
    result.peakRingStd     = peakRing.std_nm;
    result.peakRingMax     = peakRing.max_nm;
    result.peakRingMin     = peakRing.min_nm;
    result.peakRingN_px    = peakRing.n_px;
    result.peakRingMask    = peakRing.mask;
    result.Iret          = Iret;
    result.success       = true;
end


function seed = make_filled_seed(edgeMask, closeRadius)
    seed = [];
    edgeMask = bwareaopen(logical(edgeMask), 8);
    edgeMask = imclearborder(edgeMask);
    if nnz(edgeMask) < 20
        return
    end
    if closeRadius > 0
        edgeMask = imclose(edgeMask, strel('disk', closeRadius, 0));
    end
    candidate = imfill(edgeMask, 'holes');
    candidate = bwareafilt(candidate, 1);
    fraction = nnz(candidate) / numel(candidate);
    if fraction < 0.005 || fraction > 0.95
        return
    end
    seed = candidate;
end


function peakRing = compute_peak_ring_stats(Iret, BW, opts)
    depthStep_um = opts.profileDepthStep_um;
    depthAxis_um = 0:depthStep_um:opts.profileMaxDepth_um;
    nDepth = numel(depthAxis_um);

    peakRing.depthAxis_um   = depthAxis_um;
    peakRing.distProfile_nm = nan(1, nDepth);
    peakRing.depth_um       = NaN;
    peakRing.values_nm      = [];
    peakRing.mean_nm        = NaN;
    peakRing.median_nm      = NaN;
    peakRing.std_nm         = NaN;
    peakRing.max_nm         = NaN;
    peakRing.min_nm         = NaN;
    peakRing.n_px           = 0;
    peakRing.mask           = false(size(BW));

    if ~any(BW(:)) || depthStep_um <= 0 || opts.profileMaxDepth_um <= 0
        return;
    end

    per = bwperim(BW);
    D_um = bwdist(per) * opts.um_per_px;
    D_um(~BW) = NaN;

    depthBinEdges = [depthAxis_um - depthStep_um/2, depthAxis_um(end) + depthStep_um/2];
    depthBins = discretize(D_um(:), depthBinEdges);
    retVals = Iret(:);
    validD = ~isnan(depthBins) & isfinite(retVals);

    if ~any(validD)
        return;
    end

    peakRing.distProfile_nm = accumarray(depthBins(validD), retVals(validD), ...
        [nDepth 1], @local_nanmean, NaN)';

    searchMask = depthAxis_um >= opts.peakSearchDepthRange_um(1) & ...
                 depthAxis_um <= opts.peakSearchDepthRange_um(2) & ...
                 isfinite(peakRing.distProfile_nm);

    if ~any(searchMask)
        return;
    end

    searchIdx = find(searchMask);
    [~, peakRelIdx] = max(peakRing.distProfile_nm(searchIdx));
    peakIdx = searchIdx(peakRelIdx);
    peakRing.depth_um = depthAxis_um(peakIdx);

    peakPixelMask = validD & depthBins == peakIdx;
    peakVals = retVals(peakPixelMask);
    peakVals = peakVals(isfinite(peakVals));

    if isempty(peakVals)
        return;
    end

    peakRing.values_nm = peakVals;
    peakRing.mean_nm   = mean(peakVals, 'omitnan');
    peakRing.median_nm = median(peakVals, 'omitnan');
    peakRing.std_nm    = std(peakVals, 'omitnan');
    peakRing.max_nm    = max(peakVals);
    peakRing.min_nm    = min(peakVals);
    peakRing.n_px      = numel(peakVals);
    peakRing.mask      = reshape(peakPixelMask, size(BW));
end


function m = local_nanmean(x)
    m = mean(x, 'omitnan');
end


function BW = keep_largest_component(BW)
    L = bwlabel(BW, 8);
    if max(L(:)) < 1
        return;
    end
    S = regionprops(L, 'Area');
    [~, iMax] = max([S.Area]);
    BW = (L == iMax);
end


function result = make_empty_result(Iret, BW, gradMask, edgeMask, boundarySupport)
    if nargin < 3; gradMask = []; end
    if nargin < 4; edgeMask = []; end
    if nargin < 5; boundarySupport = []; end

    result.contourValues = [];
    result.contourMean   = NaN;
    result.contourStd    = NaN;
    result.contourMax    = NaN;
    result.contourMin    = NaN;
    result.xb            = [];
    result.yb            = [];
    result.xc            = NaN;
    result.yc            = NaN;
    result.R_fit         = NaN;
    result.BW            = BW;
    result.BW_initial    = BW;
    result.ridgeResponse = [];
    result.gradMask      = gradMask;
    result.edgeMask      = edgeMask;
    result.boundarySupport = boundarySupport;
    result.depthAxis_um    = [];
    result.distProfile_nm  = [];
    result.peakRingDepth_um = NaN;
    result.peakRingValues  = [];
    result.peakRingMean    = NaN;
    result.peakRingMedian  = NaN;
    result.peakRingStd     = NaN;
    result.peakRingMax     = NaN;
    result.peakRingMin     = NaN;
    result.peakRingN_px    = 0;
    result.peakRingMask    = false(size(BW));
    result.Iret          = Iret;
    result.success       = false;
end

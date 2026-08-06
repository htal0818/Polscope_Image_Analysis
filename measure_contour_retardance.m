function result = measure_contour_retardance(Iraw, opts)
% MEASURE_CONTOUR_RETARDANCE  Segment an oocyte and measure retardance at its contour.
%
%   result = measure_contour_retardance(Iraw, opts)
%
%   Takes a raw 16-bit PolScope image (double) and returns retardance values
%   sampled along the detected oocyte boundary.
%
%   INPUTS
%     Iraw  — [H x W] double, raw pixel values (e.g. 0..65535 for 16-bit)
%             Used for retardance conversion and measurement.
%     opts  — struct with fields (all optional, defaults shown):
%               retardance_ceiling_nm  (50)     Polscope ceiling in nm
%               bit_depth              (16)     image bit depth
%               sigmaBlur              (20)     Gaussian blur sigma for segmentation
%               closeRadius            (25)     morphological close disk radius
%               minArea                (5000)   minimum object area in px^2
%               boundaryInset_px       (10)     shift boundary inward onto cortex
%               thresholdMode          ('otsu') 'otsu', 'fixed', 'percentile',
%                                               or 'adaptive'
%               fixedThreshold         (500)    for 'fixed' mode
%               percentileThreshold    (30)     for 'percentile' mode
%               adaptiveSensitivity    (0.5)    for 'adaptive' mode (0-1,
%                                               higher = more foreground)
%               useGradientThreshold   (false)  add high-gradient support mask
%               gradSigma              (1.5)    Gaussian sigma for gradient image
%               gradPercentile         (92)     gradient percentile threshold
%               useEdgeThreshold       (false)  add Canny edge support mask
%               cannyThresholds        ([0.08 0.25]) Canny thresholds on normalized image
%               useBoundarySupportMask (false)  build mask from gradient/edge support
%               boundaryCloseRadius    (25)     close radius for edge/gradient mask
%               boundaryDilateRadius   (3)      dilation radius for edge/gradient support
%               useActiveContour       (false)  refine BW after thresholding
%               activeContourIterations (200)   active contour iterations
%               activeContourMethod    ('edge') 'edge' or 'Chan-Vese'
%               um_per_px              (1)      spatial calibration
%               profileMaxDepth_um     (50)     max inward depth for radial profile
%               profileDepthStep_um    (1)      inward radial bin size
%               peakSearchDepthRange_um ([0 50]) search range for peak ring
%               prevBW                 ([])     previous frame mask for fallback
%               Iseg                   ([])     separate image for segmentation
%                                               (e.g. avg of State1-4). If empty,
%                                               Iraw is used for segmentation.
%
%   OUTPUT
%     result — struct with fields:
%               contourValues  — [N x 1] retardance (nm) at each boundary point
%               contourMean    — scalar, mean contour retardance (nm)
%               contourStd     — scalar, std of contour retardance (nm)
%               contourMax     — scalar
%               contourMin     — scalar
%               xb, yb         — boundary point coordinates (after inset)
%               xc, yc         — circle-fit center
%               R_fit          — circle-fit radius (px)
%               BW             — final binary mask of the oocyte
%               BW_initial     — pre-active-contour binary mask
%               gradMask       — high-gradient support pixels
%               edgeMask       — Canny edge support pixels
%               peakRingValues — retardance values at peak radial depth
%               peakRingDepth_um — depth where outside-in profile is maximal
%               Iret           — retardance image in nm
%               success        — logical, true if boundary was found
%
%   Requires: segment_oocyte.m + circfit.m (in this repo), Image Processing Toolbox
%
%   NOTE: segmentation is delegated to segment_oocyte.m (one shared pipeline).
%   The gradient/edge/boundary-support options listed below are DEPRECATED and
%   ignored; the mask now comes from segment_oocyte with automatic egg polarity,
%   and the cortical inset follows the true inward normal (no circular assumption).

    %% Defaults
    if nargin < 2; opts = struct(); end
    def = struct( ...
        'retardance_ceiling_nm', 50, ...
        'bit_depth',             16, ...
        'sigmaBlur',             1, ...
        'closeRadius',           1, ...
        'minArea',               5000, ...
        'boundaryInset_px',      1, ...
        'thresholdMode',         'otsu', ...
        'fixedThreshold',        500, ...
        'percentileThreshold',   30, ...
        'adaptiveSensitivity',   0.5, ...
        'useGradientThreshold',  true, ...
        'gradSigma',             1.5, ...
        'gradPercentile',        92, ...
        'useEdgeThreshold',      true, ...
        'cannyThresholds',       [0.08 0.25], ...
        'edgeDilateRadius',      1, ...
        'edgeMethod',            'Canny', ...
        'gradientPercentile',    70, ...
        'useBoundarySupportMask', true, ...
        'boundaryCloseRadius',   25, ...
        'boundaryDilateRadius',  3, ...
        'useActiveContour',      false, ...
        'activeContourIterations', 200, ...
        'activeContourMethod',   'edge', ...
        'activeSmoothFactor',    1.0, ...
        'activeContractionBias', 0.0, ...
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

    %% Boundary detection — delegated to the shared segment_oocyte pipeline
    % Segment on a separate image if provided (e.g. the State1..4 sum, which has
    % far better egg/background contrast than the retardance image), then copy
    % that mask onto the retardance image for measurement. If no Iseg is given,
    % segment the retardance image itself. All segmentation lives in ONE place
    % (segment_oocyte.m), shared with the flow scripts.
    if ~isempty(opts.Iseg)
        Iseg = opts.Iseg;
    else
        Iseg = Iraw;
    end

    segParams = struct( ...
        'sigmaBlur',               opts.sigmaBlur, ...
        'closeRadius',             opts.closeRadius, ...
        'minArea',                 opts.minArea, ...
        'thresholdMode',           opts.thresholdMode, ...
        'fixedThreshold',          opts.fixedThreshold, ...
        'percentileThreshold',     opts.percentileThreshold, ...
        'adaptSensitivity',        opts.adaptiveSensitivity, ...
        'edgeMethod',              opts.edgeMethod, ...
        'edgeDilateRadius',        opts.edgeDilateRadius, ...
        'gradientPercentile',      opts.gradientPercentile, ...
        'useActiveContour',        opts.useActiveContour, ...
        'activeContourIterations', opts.activeContourIterations, ...
        'activeContourMethod',     opts.activeContourMethod, ...
        'useCaching',              false, ...
        'smoothMask',              false, ...
        'smoothBoundary',          false);

    % Egg polarity:
    %   - four-state solid disc (Iseg given): leave segFromMask UNSET so
    %     segment_oocyte auto-detects whether the egg is darker or brighter.
    %   - retardance image (no Iseg): the egg is defined by its bright cortical
    %     rim, so force the bright class (segFromMask = false); segment_oocyte
    %     then fills the disc enclosed by the rim.
    %   - a caller may override either by setting opts.segFromMask explicitly.
    if isfield(opts, 'segFromMask') && ~isempty(opts.segFromMask)
        segParams.segFromMask = opts.segFromMask;
    elseif isempty(opts.Iseg)
        segParams.segFromMask = false;
    end

    [BW, xc, yc, R_fit, polyXY] = segment_oocyte(Iseg, segParams, []);

    % Legacy support masks: no longer used for detection, kept (empty) so the
    % result struct and the overlay saver stay backward-compatible.
    gradMask = false(H, W);
    edgeMask = false(H, W);
    boundarySupport = false(H, W);

    if isempty(polyXY) && ~isempty(opts.prevBW)
        BW = opts.prevBW;
    end
    if ~any(BW(:))
        result = make_empty_result(Iret, BW, gradMask, edgeMask, boundarySupport);
        return;
    end
    BW_initial = BW;

    %% Boundary contour (the TRUE outline from the mask — NOT assumed circular)
    if ~isempty(polyXY)
        xb = polyXY(:,1);
        yb = polyXY(:,2);
    else
        B = bwboundaries(BW);
        if isempty(B)
            result = make_empty_result(Iret, BW, gradMask, edgeMask, boundarySupport);
            return;
        end
        [~, iLongest] = max(cellfun(@(p) size(p,1), B));
        bnd = B{iLongest};
        yb = bnd(:,1);
        xb = bnd(:,2);
    end

    %% Inset the boundary inward along the TRUE inward normal (shape-agnostic)
    % Peel each boundary point in by boundaryInset_px along the local inward
    % normal, taken from the gradient of the in-mask distance transform. This
    % follows dents and elongation of a non-circular egg, unlike the old method
    % that shrank every point radially toward a single circle-fit centre.
    Din = bwdist(~BW);                       % depth inside the egg (px), 0 at edge
    [Gx, Gy] = imgradientxy(Din, 'central'); % Gx = x/column gradient, Gy = y/row
    F_Gx = griddedInterpolant({1:H, 1:W}, Gx, 'linear', 'nearest');
    F_Gy = griddedInterpolant({1:H, 1:W}, Gy, 'linear', 'nearest');
    nx = F_Gx(yb, xb);                       % x-component of inward normal
    ny = F_Gy(yb, xb);                       % y-component of inward normal
    nmag = hypot(nx, ny) + eps;
    nx = nx ./ nmag;                         % inward = +gradient of interior depth
    ny = ny ./ nmag;

    % Guarantee the normal points into the egg (toward the centroid), even where
    % the distance-transform gradient is degenerate.
    S = regionprops(BW, 'Centroid');
    if ~isempty(S)
        cen = S(1).Centroid;
        flip = ((cen(1) - xb) .* nx + (cen(2) - yb) .* ny) < 0;
        nx(flip) = -nx(flip);
        ny(flip) = -ny(flip);
    end

    xb = xb + opts.boundaryInset_px * nx;
    yb = yb + opts.boundaryInset_px * ny;

    %% Interpolate retardance (nm) at the inset boundary locations
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

function BW_boundary = make_boundary_support_mask(boundarySupport, opts)
    BW_boundary = boundarySupport;

    if opts.edgeDilateRadius > 0
        BW_boundary = imdilate(BW_boundary, strel('disk', opts.edgeDilateRadius));
    end

    if opts.boundaryDilateRadius > 0
        BW_boundary = imdilate(BW_boundary, strel('disk', opts.boundaryDilateRadius));
    end

    BW_boundary = imclose(BW_boundary, strel('disk', opts.boundaryCloseRadius));
    BW_boundary = imfill(BW_boundary, 'holes');
    BW_boundary = bwareaopen(BW_boundary, opts.minArea);
    BW_boundary = keep_largest_component(BW_boundary);
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
    result.gradMask      = gradMask;
    result.edgeMask      = edgeMask;
    result.boundarySupport = boundarySupport;
    result.Iret          = Iret;
    result.success       = false;
end

function I_norm = normalize01(I)
    I = double(I);
    lo = min(I(:));
    hi = max(I(:));
    if hi > lo
        I_norm = (I - lo) ./ (hi - lo);
    else
        I_norm = zeros(size(I));
    end
end

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
%               BW             — binary mask of the oocyte
%               Iret           — retardance image in nm
%               success        — logical, true if boundary was found
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
        'boundaryInset_px',      10, ...
        'thresholdMode',         'otsu', ...
        'fixedThreshold',        500, ...
        'percentileThreshold',   30, ...
        'prevBW',                [], ...
        'Iseg',                  [], ...
        'adaptiveSensitivity',   0.5);
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

    %% Boundary detection
    % Use separate segmentation image if provided (e.g. sum of State1-4),
    % otherwise fall back to Iraw.
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

    % Fallback: gradient-based
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
        % Failed — return empty result
        result = make_empty_result(Iret, BW);
        return;
    end

    %% Extract boundary contour
    B = bwboundaries(BW);
    if isempty(B)
        result = make_empty_result(Iret, BW);
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
    result.Iret          = Iret;
    result.success       = true;
end

function result = make_empty_result(Iret, BW)
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
    result.Iret          = Iret;
    result.success       = false;
end

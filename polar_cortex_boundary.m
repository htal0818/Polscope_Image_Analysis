function [R_theta, xc, yc, info] = polar_cortex_boundary(Iret, BW_thresh, params)
% POLAR_CORTEX_BOUNDARY  Per-frame cortex boundary by radial peak search.
%
% Each frame's boundary is estimated independently — no temporal state,
% no per-frame initial condition, no drift. Suited to oocyte cortex
% imaging where the cortex is the dominant outer radial feature.
%
% Algorithm:
%   1. Centroid (xc, yc) from BW_thresh (area-weighted).
%   2. Coarse cortex radius R0 = sqrt(area / pi).
%   3. Polar sample grid: nTheta angles, nR radii in
%      [cortexSearchMinFrac, cortexSearchMaxFrac] * R0.
%   4. Sample Iret_smoothed in polar coords via griddedInterpolant.
%   5. Per angle: findpeaks with prominence/height filters, keep the
%      strongest local maximum. Fall back to global max if no peaks.
%   6. Circular median (or chosen smoothdata method) on R(theta).
%   7. Wrap-around fill of any NaN angles via interp1.
%
% Inputs:
%   Iret      - [H x W] retardance image (nm)
%   BW_thresh - [H x W] logical, coarse oocyte mask (provides center)
%   params    - struct with fields:
%       nTheta              - number of angular samples (e.g. 720)
%       nR                  - number of radial samples (e.g. 400)
%       cortexSearchMinFrac - inner search bound (fraction of R0)
%       cortexSearchMaxFrac - outer search bound (fraction of R0)
%       smoothSigma         - Gaussian sigma (px) on Iret before sampling
%       minPeakValue        - findpeaks 'MinPeakHeight'
%       minPeakProminence   - findpeaks 'MinPeakProminence'
%       smoothMethod        - smoothdata kernel ('movmedian', 'sgolay', ...)
%       smoothWindow        - smoothdata window size (angular samples)
%
% Outputs:
%   R_theta - [1 x nTheta] cortex radius per angle (NaN if not found)
%   xc, yc  - centroid used (image coords, 1-based)
%   info    - struct with diagnostic fields (nPeaksFound, nFallback, etc.)

info = struct('nPeaksFound', 0, 'nFallback', 0, 'nMissing', 0, ...
              'nContinuityRevised', 0, ...
              'R0_px', NaN, 'Rmin_px', NaN, 'Rmax_px', NaN);

R_theta = [];

S = regionprops(BW_thresh, 'Centroid', 'Area');
if isempty(S)
    return;
end
[~, iMax] = max([S.Area]);
xc = S(iMax).Centroid(1);
yc = S(iMax).Centroid(2);

R0 = sqrt(S(iMax).Area / pi);
Rmin_px = max(1, params.cortexSearchMinFrac * R0);
Rmax_px = params.cortexSearchMaxFrac * R0;
info.R0_px   = R0;
info.Rmin_px = Rmin_px;
info.Rmax_px = Rmax_px;

[H, W] = size(Iret);

theta = linspace(0, 2*pi, params.nTheta + 1);
theta(end) = [];
r = linspace(Rmin_px, Rmax_px, params.nR);

% Light denoise on the measurement image, then sample in polar coords.
Iret_s = imgaussfilt(Iret, params.smoothSigma);
F = griddedInterpolant({1:H, 1:W}, Iret_s, 'linear', 'nearest');

[Rgrid, Tgrid] = ndgrid(r, theta);
Xq = xc + Rgrid .* cos(Tgrid);
Yq = yc + Rgrid .* sin(Tgrid);

inBounds = Xq >= 1 & Xq <= W & Yq >= 1 & Yq <= H;
Ipolar = nan(size(Xq));
Ipolar(inBounds) = F(Yq(inBounds), Xq(inBounds));

% Pass 1: collect ALL candidate peaks per ray, don't pick yet.
nT = params.nTheta;
candRadii = cell(1, nT);
candVals  = cell(1, nT);

for j = 1:nT
    profile = Ipolar(:, j);
    if all(isnan(profile)) || max(profile, [], 'omitnan') <= 0
        info.nMissing = info.nMissing + 1;
        candRadii{j} = [];
        candVals{j}  = [];
        continue;
    end

    profile_clean = profile;
    profile_clean(isnan(profile_clean)) = 0;
    try
        [pkVals, pkLocs] = findpeaks(profile_clean, ...
            'MinPeakHeight',     params.minPeakValue, ...
            'MinPeakProminence', params.minPeakProminence);
    catch
        pkVals = []; pkLocs = [];
    end

    if isempty(pkLocs)
        % Fallback: global max becomes the single candidate.
        [pkVal, pkLoc] = max(profile_clean);
        if pkVal >= params.minPeakValue
            candRadii{j} = r(pkLoc);
            candVals{j}  = pkVal;
            info.nFallback = info.nFallback + 1;
        else
            candRadii{j} = [];
            candVals{j}  = [];
            info.nMissing = info.nMissing + 1;
        end
    else
        candRadii{j} = r(pkLocs);
        candVals{j}  = pkVals;
        info.nPeaksFound = info.nPeaksFound + 1;
    end
end

% Pass 2: initial pick = outermost significant peak per ray. Used to
% seed the continuity refinement; no NaN bridging.
R_theta = nan(1, nT);
for j = 1:nT
    cV = candVals{j};
    cR = candRadii{j};
    if isempty(cR); continue; end
    valid = cV > params.peakKeepFrac * max(cV);
    if any(valid)
        R_theta(j) = cR(find(valid, 1, 'last'));
    else
        [~, iMax] = max(cV);
        R_theta(j) = cR(iMax);
    end
end

% Pass 3: angular-continuity refinement. For each ray, recompute the
% pick restricted to candidates within maxJumpPx of the local circular
% median. Within that band, still prefer the outermost significant
% peak so polar-body bulges are kept. No NaN, no interpolation: every
% ray is assigned the radius of an actual candidate peak that exists
% in its own profile and is consistent with its neighborhood.
info.nContinuityRevised = 0;
if isfield(params, 'useAngularContinuity') && params.useAngularContinuity
    R_pad = [R_theta R_theta R_theta];
    R_med = smoothdata(R_pad, 'movmedian', params.continuityMedianWin, ...
                       'includenan');
    R_med = R_med(nT + 1 : 2*nT);
    R_med_global = median(R_theta(~isnan(R_theta)), 'omitnan');

    for j = 1:nT
        cV = candVals{j};
        cR = candRadii{j};
        if isempty(cR); continue; end

        R_pred = R_med(j);
        if isnan(R_pred); R_pred = R_med_global; end
        if isnan(R_pred); continue; end

        within = abs(cR - R_pred) <= params.maxJumpPx;
        if ~any(within); continue; end   % keep pass-2 pick

        cR_w = cR(within);
        cV_w = cV(within);
        valid = cV_w > params.peakKeepFrac * max(cV_w);
        if any(valid)
            newR = cR_w(find(valid, 1, 'last'));
        else
            [~, iMax] = max(cV_w);
            newR = cR_w(iMax);
        end

        if ~isnan(R_theta(j)) && abs(newR - R_theta(j)) > 0
            info.nContinuityRevised = info.nContinuityRevised + 1;
        end
        R_theta(j) = newR;
    end
end

% Periodic smoothing — pad with copies, smooth, trim.
R_ext = [R_theta R_theta R_theta];
R_ext = smoothdata(R_ext, params.smoothMethod, params.smoothWindow, ...
                   'includenan');
R_theta = R_ext(nT + 1 : 2*nT);

% Fill any remaining NaN via circular interpolation.
nanIdx = isnan(R_theta);
if any(nanIdx) && any(~nanIdx)
    valIdx  = find(~nanIdx);
    nanLocs = find(nanIdx);
    x_ext = [valIdx - nT, valIdx, valIdx + nT];
    y_ext = [R_theta(valIdx), R_theta(valIdx), R_theta(valIdx)];
    R_theta(nanIdx) = interp1(x_ext, y_ext, nanLocs, 'linear', 'extrap');
end

% Final low-pass smoothing — Savitzky-Golay preserves the polar body
% bulge (it's a localized polynomial fit, not a kernel average) while
% removing per-angle sawtooth oscillations from peak-pick noise.
if isfield(params, 'sgolayWindow') && params.sgolayWindow > 0
    R_ext = [R_theta R_theta R_theta];
    R_ext = smoothdata(R_ext, 'sgolay', params.sgolayWindow);
    R_theta = R_ext(nT + 1 : 2*nT);
end

end

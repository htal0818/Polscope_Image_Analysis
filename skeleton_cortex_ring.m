function [R_theta, xc, yc, info] = skeleton_cortex_ring(Iret, BW_thresh, params)
% SKELETON_CORTEX_RING  Find the cortex boundary by skeletonizing the bright
% retardance ring around the oocyte.
%
% Fundamentally different from peak-picking: the cortex is identified as
% a CONNECTED geometric structure in the image -- the largest bright ring
% within a band of the rough oocyte mask. External bright objects (debris,
% polar bodies, separate cells) are disconnected from the cortex ring by
% dark medium and get rejected automatically by bwareafilt.
%
% Algorithm:
%   1. Find bright pixels likely to be cortex: Iret > prctile(Iret(BW_thresh), pctLevel)
%   2. Restrict to a band around the rough oocyte boundary (within bandWidth_um
%      of BW_thresh's edge). Excludes interior bright noise.
%   3. Keep only the largest connected component (= the cortex ring).
%   4. Skeletonize to get the ring centerline.
%   5. Convert skeleton (x,y) -> (theta, r) about the centroid.
%   6. Bin into nTheta angles, take median radius per bin.
%   7. Fill missing angles via circular interp1; light circular smoothing.
%
% Inputs:
%   Iret      - [H x W] retardance image (nm or raw pixel units)
%   BW_thresh - [H x W] logical, rough oocyte mask (provides centroid + band)
%   params    - struct with fields:
%       nTheta              - number of angular samples (e.g. 720)
%       pctLevel            - percentile of Iret(BW_thresh) for cortex threshold (e.g. 90)
%       bandWidth_um        - microns wide band around BW_thresh edge to search
%       px_per_um           - microns -> pixels conversion
%       smoothWindow        - circular median smoothing window (angular samples)
%
% Outputs:
%   R_theta - [1 x nTheta] cortex radius per angle (NaN-free)
%   xc, yc  - centroid used
%   info    - struct with diagnostic fields

info = struct('cortexLevel', NaN, 'nRingPixels', 0, 'nSkeletonPixels', 0, ...
              'nMissingAngles', 0);

R_theta = [];

S = regionprops(BW_thresh, 'Centroid', 'Area');
if isempty(S)
    return;
end
[~, iMax] = max([S.Area]);
xc = S(iMax).Centroid(1);
yc = S(iMax).Centroid(2);

[H, W] = size(Iret);

% Step 1: cortex-level threshold on Iret
%   Use only pixels inside the rough oocyte mask to estimate the level so
%   the medium / external objects don't influence the percentile.
cortexLevel = prctile(Iret(BW_thresh), params.pctLevel);
brightPixels = Iret > cortexLevel;
info.cortexLevel = cortexLevel;

% Step 2: restrict to band around BW_thresh edge
%   D = signed distance: positive inside, negative outside
D_in  = bwdist(~BW_thresh);     % distance to outside (positive inside)
D_out = bwdist(BW_thresh);      % distance to inside  (positive outside)
band_px = round(params.bandWidth_um * params.px_per_um);
nearBoundary = (D_in > 0 & D_in < band_px) | (D_out > 0 & D_out < band_px);
ringMask = brightPixels & nearBoundary;

% Step 3: keep largest connected component (the cortex ring)
if ~any(ringMask(:))
    return;
end
ringMask = bwareafilt(ringMask, 1);
info.nRingPixels = nnz(ringMask);

% Step 4: skeletonize
skel = bwskel(ringMask);
info.nSkeletonPixels = nnz(skel);

if info.nSkeletonPixels < 10
    return;
end

% Step 5: (x, y) -> (theta, r)
[ys, xs] = find(skel);
theta_pts = atan2(ys - yc, xs - xc);
theta_pts(theta_pts < 0) = theta_pts(theta_pts < 0) + 2*pi;
r_pts = hypot(xs - xc, ys - yc);

% Step 6: bin by angle, median radius per bin
nT = params.nTheta;
binEdges = linspace(0, 2*pi, nT + 1);
[~, ~, binIdx] = histcounts(theta_pts, binEdges);

R_theta = nan(1, nT);
for j = 1:nT
    inBin = binIdx == j;
    if any(inBin)
        R_theta(j) = median(r_pts(inBin));
    end
end
info.nMissingAngles = nnz(isnan(R_theta));

% Step 7a: fill missing angles via circular interp1
nanIdx = isnan(R_theta);
if any(nanIdx) && any(~nanIdx)
    valIdx  = find(~nanIdx);
    nanLocs = find(nanIdx);
    x_ext = [valIdx - nT, valIdx, valIdx + nT];
    y_ext = [R_theta(valIdx), R_theta(valIdx), R_theta(valIdx)];
    R_theta(nanIdx) = interp1(x_ext, y_ext, nanLocs, 'linear', 'extrap');
end

% Step 7b: light circular median smoothing
if isfield(params, 'smoothWindow') && params.smoothWindow > 1
    R_pad = [R_theta R_theta R_theta];
    R_pad = smoothdata(R_pad, 'movmedian', params.smoothWindow);
    R_theta = R_pad(nT + 1 : 2*nT);
end

end

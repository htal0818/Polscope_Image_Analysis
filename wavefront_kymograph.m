% wavefront_kymograph.m
% Per-angle (depth x time) carpet view of inward retardance.
%
% Why this is a separate script from contour_retardance.m:
%   contour_retardance.m's distKymo (depth x time) angle-averages every
%   interior pixel at each depth, which washes out angle-localized internal
%   signals: a feature in 1/N of the cortex contributes ~1/N to the average
%   while the bright cortex band dominates everywhere. The user-confirmed
%   wash-out shows only the cortex stripe; known waves at ~25 um depth
%   never appear in distKymo.
%
%   This script builds the full (angle x depth x frame) cube once, then
%   stacks every angle's (depth x time) kymograph vertically into a single
%   "carpet" image. No reductions, no angle-picking, no thresholding —
%   the user inspects the carpet directly to see which angles carry the
%   wave and how it propagates around the cortex.
%
% Carpet structure:
%   - Y axis: (angle, depth) flattened. Each contiguous block of nDepth
%     rows is one angle's depth-vs-time panel, with depth=0 (cortex) at
%     the block's top row and increasing depth going down within the
%     block. Angle blocks are stacked from theta=0 (top) to theta=2*pi
%     (bottom).
%   - X axis: time (shared across all angles).
%   - Pixel value: retardance at that (angle, depth, time) in nm.
%   - Faint horizontal separators every nDepth rows mark angle blocks.
%   - Y-axis ticks at block centres are labelled with theta in degrees
%     (every Nth angle to avoid crowding).
%
% A wave in one cortex sector appears as a localized streak below one
% (or a few contiguous) cortex-stripe rows. With nAngleBins=100 and
% ~10 underlying cortex sectors the wave shows up as a "thick" band of
% ~10 contiguous similar streaks — redundancy that aids detection
% rather than diluting it (each sector's panel is independently kept).
%
% The per-frame (angle x depth) map is computed in one vectorized
% accumarray pass over the interior pixels using the distance transform.
% Inputs are the same image stacks consumed by contour_retardance.m
% (4-state Polscope is the default; retardance and multipage are also
% supported).
%
% Outputs (in outDir):
%   carpet_per_angle.png/.fig          per-angle (depth x time) carpet
%   wavefront_kymograph_results.mat    full cube + carpet + axes + diagnostics
%
% REQUIREMENTS:
%   - Image Processing Toolbox
%   - circfit.m (in this repo)

clear all; close all; clc

%% ========================== USER INPUTS ==================================
% --- Input mode (mirrors contour_retardance.m) ---
%   'retardance' : folder of pre-computed retardance TIFFs
%   'four_state' : raw 4-state Polscope (State1..State4) — default
%   'multipage'  : single multipage TIFF
inputMode = 'four_state';

base_dir = '/Users/hridaytalreja/Desktop/Mar_2026_data/2026_04_01/SMS_2026_0401_1040_1/Pos0/';
retardance_pattern = '*1_Retardance*';

% Mask source pattern (e.g. '*State1*'); '' to segment from retardance.
% Ignored in four_state mode (the State1..State4 sum is used automatically).
mask_pattern = 'four_state';

% multipage_path = '/path/to/stack.tif';

% --- Optional crop ---
doCrop = false;
cropRect = [500 500 2000 2000];

% --- Calibration ---
dt_sec    = 15;
px_per_um = 6.25/2;
retardance_ceiling_nm = 50;
bit_depth = 16;

% --- Segmentation: edge-detection + morphology on the 4-state sum.
%     Defaults match contour_retardance.m. ---
thresholdMode    = 'edge';     % 'edge' or 'otsu'
edgeMethod       = 'Canny';    % edge() method: 'Canny', 'Sobel', 'Prewitt', 'log'
edgeDilateRadius = 1;          % dilation radius (px) to close edge gaps before fill
sigmaBlur        = 1;          % Gaussian blur sigma (px)
closeRadius      = 1;          % morphological close disk radius (px)
minArea          = 5000;       % minimum component area (px^2)

% --- Angle x depth grid ---
nAngleBins   = 100;     % angular bins around centroid (covers 0..2*pi)
depthStep_um = 0.5;     % depth bin step (microns)
maxDepth_um  = 50;      % how far inward from cortex to sample (microns)

% --- Carpet plot appearance ---
carpetAngleLabelEvery = 10;     % label every Nth angle on the Y-axis
carpetSeparatorAlpha  = 0.25;   % opacity of the inter-angle-block separator lines
carpetClimQuantiles   = [0.02 0.98];  % robust color limits

% --- Output ---
outDir = fullfile(fileparts(mfilename('fullpath')), 'wavefront_kymograph_out');

%% ========================== LOAD IMAGE LIST ===============================

switch inputMode
    case 'retardance'
        d = dir(fullfile(base_dir, retardance_pattern));
        if isempty(d)
            error('No retardance images found matching "%s" in %s', ...
                retardance_pattern, base_dir);
        end
        [~, sortIdx] = sort({d.name});
        d = d(sortIdx);
        nFrames = numel(d);
        readFrame = @(t) double(imread(fullfile(d(t).folder, d(t).name)));

    case 'four_state'
        state_patterns = {'*State1*', '*State2*', '*State3*', '*State4*'};
        ds = cell(1,4);
        for si = 1:4
            ds{si} = dir(fullfile(base_dir, state_patterns{si}));
            [~, idx] = sort({ds{si}.name});
            ds{si} = ds{si}(idx);
        end
        readMask = @(t) double(imread(fullfile(ds{1}(t).folder, ds{1}(t).name))) ...
                      + double(imread(fullfile(ds{2}(t).folder, ds{2}(t).name))) ...
                      + double(imread(fullfile(ds{3}(t).folder, ds{3}(t).name))) ...
                      + double(imread(fullfile(ds{4}(t).folder, ds{4}(t).name)))/4;

        d = dir(fullfile(base_dir, retardance_pattern));
        if isempty(d)
            error('No retardance images found matching "%s" in %s', ...
                retardance_pattern, base_dir);
        end
        [~, sortIdx] = sort({d.name});
        d = d(sortIdx);

        nFrames     = min([numel(d), cellfun(@numel, ds)]);
        nMaskFrames = min(cellfun(@numel, ds));
        readFrame   = @(t) double(imread(fullfile(d(t).folder, d(t).name)));
        useMaskSource = true;

    case 'multipage'
        info = imfinfo(multipage_path);
        nFrames = numel(info);
        readFrame = @(t) double(imread(multipage_path, t));

    otherwise
        error('Unknown inputMode: %s', inputMode);
end

if ~exist('useMaskSource', 'var'); useMaskSource = false; end
if ~useMaskSource && ~isempty(mask_pattern)
    ds_mask = dir(fullfile(base_dir, mask_pattern));
    if ~isempty(ds_mask)
        [~, mIdx] = sort({ds_mask.name});
        ds_mask = ds_mask(mIdx);
        nMaskFrames = numel(ds_mask);
        readMask = @(t) double(imread(fullfile(ds_mask(t).folder, ds_mask(t).name)));
        useMaskSource = true;
    end
end

fprintf('Found %d frames (mode: %s)\n', nFrames, inputMode);

%% ========================== SETUP =========================================
if ~exist(outDir, 'dir')
    [status, msg] = mkdir(outDir);
    if ~status
        error('Could not create output directory "%s": %s', outDir, msg);
    end
end

um_per_px = 1 / px_per_um;
maxPixVal = 2^bit_depth - 1;

time_sec = (0:nFrames-1)' * dt_sec;
time_min = time_sec / 60;

% Angle bins: 0..2*pi
angleEdges   = linspace(0, 2*pi, nAngleBins+1);
angleCenters = (angleEdges(1:end-1) + angleEdges(2:end)) / 2;

% Depth bins: 0 (cortex) .. maxDepth_um
depthEdges  = 0 : depthStep_um : maxDepth_um;
depthAxis_um = (depthEdges(1:end-1) + depthEdges(2:end)) / 2;
nDepth = numel(depthAxis_um);

% Full (angle x depth x frame) cube — every per-frame slice is stored intact.
sliceCube    = nan(nAngleBins, nDepth, nFrames);
centroidXY   = nan(nFrames, 2);
fitRadius_um = nan(nFrames, 1);

cache_prevBW = [];

%% ========================== MAIN LOOP =====================================
fprintf('Computing per-frame inward profiles...\n');
tic;

for fr = 1:nFrames
    Iraw = readFrame(fr);
    if doCrop; Iraw = imcrop(Iraw, cropRect); end
    [H, W] = size(Iraw); %#ok<ASGLU>

    Iret = (Iraw / maxPixVal) * retardance_ceiling_nm;

    % --- Choose segmentation source ---
    if useMaskSource && fr <= nMaskFrames
        Iseg = readMask(fr);
        if doCrop; Iseg = imcrop(Iseg, cropRect); end
        segFromMask = true;
    else
        Iseg = Iraw;
        segFromMask = false;
    end

    % --- Segmentation ---
    I_blur = imgaussfilt(Iseg, sigmaBlur);
    I_norm = I_blur / max(I_blur(:));

    switch thresholdMode
        case 'edge'
            edges = edge(I_norm, edgeMethod);
            se_edge = strel('disk', edgeDilateRadius);
            edges = imdilate(edges, se_edge);
            BW = imfill(edges, 'holes');
            % If we segmented from a "dark oocyte" 4-state sum and the
            % filled mask covers most of the frame, invert.
            if segFromMask && sum(BW(:)) > 0.5 * numel(BW)
                BW = ~BW;
            end

        case 'otsu'
            Totsu  = graythresh(I_norm);
            if segFromMask
                BW = I_norm < Totsu;
            else
                BW = I_norm > Totsu;
            end

        otherwise
            error('Unknown thresholdMode: %s', thresholdMode);
    end

    se = strel('disk', closeRadius);
    BW = imclose(BW, se);
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, minArea);

    % Keep largest component
    L = bwlabel(BW, 8);
    if max(L(:)) >= 1
        S = regionprops(L, 'Area');
        [~, iMax] = max([S.Area]);
        BW = (L == iMax);
    elseif ~isempty(cache_prevBW)
        BW = cache_prevBW;
    else
        fprintf('  Frame %d: no boundary found, skipping.\n', fr);
        continue;
    end
    cache_prevBW = BW;

    % --- Boundary -> circle fit (centroid + radius for diagnostics) ---
    Bnd = bwboundaries(BW);
    if isempty(Bnd); continue; end
    [~, iLongest] = max(cellfun(@(p) size(p,1), Bnd));
    bnd = Bnd{iLongest};
    [R_fit, xc, yc] = circfit(bnd(:,2), bnd(:,1));
    centroidXY(fr,:) = [xc, yc];
    fitRadius_um(fr) = R_fit * um_per_px;

    % --- Distance transform: depth-from-cortex for every interior pixel ---
    D = bwdist(~BW) * um_per_px;

    % Restrict to interior pixels within maxDepth_um of the cortex.
    interior = BW & (D <= maxDepth_um);
    if ~any(interior(:)); continue; end

    [yIdx, xIdx] = find(interior);
    theta = mod(atan2(yIdx - yc, xIdx - xc), 2*pi);

    aBin = discretize(theta,         angleEdges);
    dBin = discretize(D(interior),   depthEdges);

    valid = ~isnan(aBin) & ~isnan(dBin);
    if ~any(valid); continue; end

    linIdx  = sub2ind(size(interior), yIdx(valid), xIdx(valid));
    retVals = Iret(linIdx);

    % Per-frame (angle x depth) mean retardance map (single accumarray).
    slice = accumarray( ...
        [aBin(valid), dBin(valid)], retVals, ...
        [nAngleBins, nDepth], @mean, NaN);

    sliceCube(:, :, fr) = slice;

    if mod(fr, max(1, round(nFrames/10))) == 0
        fprintf('  Processed %d / %d frames (%.0f%%)\n', fr, nFrames, 100*fr/nFrames);
    end
end

elapsed = toc;
fprintf('Done! %.1f sec total (%.2f sec/frame)\n\n', elapsed, elapsed/nFrames);

%% ========================== CARPET ASSEMBLY ===============================
% Reshape (nAngleBins x nDepth x nFrames) cube into a 2D carpet image.
% After permute, the in-memory axis order is (depth, angle, frame).
% Column-major reshape over the first two dims stacks angle blocks
% vertically with depth contiguous within each block: rows 1..nDepth =
% angle 1 (depth 0 at top, increasing depth down), rows nDepth+1..2*nDepth
% = angle 2, and so on.
carpet = reshape(permute(sliceCube, [2 1 3]), ...
                 [nDepth * nAngleBins, nFrames]);   % rows x time-cols

%% ========================== PLOT ==========================================
fprintf('Generating carpet plot...\n');

fig = figure('Position', [80 60 900 1400]);
ax  = axes('Parent', fig);
imagesc(ax, time_min, 1:size(carpet,1), carpet);
set(ax, 'YDir', 'reverse');                   % depth=0 at top of each block
xlabel(ax, 'Time (min)', 'FontSize', 12);
ylabel(ax, 'Angle around cortex (deg)  \rightarrow  depth within block (\mum)', ...
       'FontSize', 12);
title(ax, 'Per-angle inward retardance — carpet view', 'FontSize', 14);

% Robust color limits over the whole carpet
finiteVals = carpet(isfinite(carpet));
if ~isempty(finiteVals)
    qLo = quantile(finiteVals, carpetClimQuantiles(1));
    qHi = quantile(finiteVals, carpetClimQuantiles(2));
    if qHi > qLo
        caxis(ax, [qLo qHi]);
    end
end
colormap(ax, parula);
cb = colorbar(ax);
cb.Label.String   = 'Retardance (nm)';
cb.Label.FontSize = 11;

% Faint horizontal separators between angle blocks
hold(ax, 'on');
for a = 1:(nAngleBins - 1)
    yLine = a * nDepth + 0.5;
    line(ax, ax.XLim, [yLine yLine], ...
         'Color', [0 0 0 carpetSeparatorAlpha], ...
         'LineWidth', 0.5);
end
hold(ax, 'off');

% Y-ticks at block centres, labelled with theta in degrees.
blockCentres = ((1:nAngleBins) - 0.5) * nDepth + 0.5;
labelMask    = mod(0:nAngleBins-1, carpetAngleLabelEvery) == 0;
ax.YTick      = blockCentres(labelMask);
ax.YTickLabel = arrayfun(@(rad) sprintf('%.0f\\circ', rad2deg(rad)), ...
                         angleCenters(labelMask), ...
                         'UniformOutput', false);

exportgraphics(fig, fullfile(outDir, 'carpet_per_angle.png'), 'Resolution', 200);
savefig(fig, fullfile(outDir, 'carpet_per_angle.fig'));
close(fig);

%% ========================== SAVE DATA =====================================
results = struct();
results.sliceCube           = sliceCube;             % full (angle x depth x frame) cube
results.carpet              = carpet;                % carpet image (rows=angle*depth, cols=time)
results.angleCenters_rad    = angleCenters;
results.angleEdges_rad      = angleEdges;
results.depthAxis_um        = depthAxis_um;
results.depthEdges_um       = depthEdges;
results.time_sec            = time_sec;
results.time_min            = time_min;
results.centroidXY          = centroidXY;
results.fitRadius_um        = fitRadius_um;
results.nFrames             = nFrames;
results.dt_sec              = dt_sec;
results.px_per_um           = px_per_um;
results.retardance_ceiling_nm = retardance_ceiling_nm;
results.bit_depth           = bit_depth;
results.nAngleBins          = nAngleBins;
results.depthStep_um        = depthStep_um;
results.maxDepth_um         = maxDepth_um;
results.thresholdMode       = thresholdMode;
results.edgeMethod          = edgeMethod;
results.inputMode           = inputMode;
results.carpetAngleLabelEvery = carpetAngleLabelEvery;
results.carpetSeparatorAlpha  = carpetSeparatorAlpha;
results.carpetClimQuantiles   = carpetClimQuantiles;

save(fullfile(outDir, 'wavefront_kymograph_results.mat'), '-struct', 'results');
fprintf('Saved results to %s\n', fullfile(outDir, 'wavefront_kymograph_results.mat'));

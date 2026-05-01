% wavefront_kymograph.m
% Per-frame max-peak inward retardance profile -> depth-vs-time kymograph.
%
% Why this is a separate script from contour_retardance.m:
%   contour_retardance.m angle-averages the inward retardance profile, which
%   washes out the localized, weak alignment band that propagates inward
%   during the contraction wave. This script keeps the angle resolution per
%   frame: it builds an (angle x depth) retardance map, picks the angular
%   direction whose inward profile peaks the highest (global peak when
%   scanning inward), and stacks that single 1D profile for each frame into
%   a kymograph.
%
% The map per frame is computed in one vectorized accumarray pass over the
% interior pixels using the distance transform — no per-angle ray loop, no
% per-boundary-point loop. Only the chosen 1D profile is kept (no 3D cube).
%
% Inputs are the same image stacks consumed by contour_retardance.m
% (4-state Polscope is the default; retardance and multipage are also
% supported).
%
% Outputs (in outDir):
%   inward_kymograph_max_peak.png/.fig    1D inward profile (chosen angle) per frame
%   inward_kymograph_avg.png/.fig         angle-averaged profile per frame (for comparison)
%   wave_origin_angle_vs_time.png/.fig    chosen direction over time
%   wave_peak_depth_vs_time.png/.fig      depth of the band over time -> wave speed
%   wavefront_kymograph_results.mat       per-frame profiles + selections
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

base_dir = '/path/to/your/data/Pos0/';
retardance_pattern = '*1_Retardance*';

% Mask source pattern (e.g. '*State1*'); '' to segment from retardance.
% Ignored in four_state mode (the State1..State4 sum is used automatically).
mask_pattern = '';

% multipage_path = '/path/to/stack.tif';

% --- Optional crop ---
doCrop = false;
cropRect = [500 500 2000 2000];

% --- Calibration ---
dt_sec    = 15;
px_per_um = 6.25;
retardance_ceiling_nm = 50;
bit_depth = 16;

% --- Segmentation: edge-detection + morphology on the 4-state sum.
%     Defaults match contour_retardance.m. ---
thresholdMode    = 'edge';     % 'edge' or 'otsu'
edgeMethod       = 'Canny';    % edge() method: 'Canny', 'Sobel', 'Prewitt', 'log'
edgeDilateRadius = 3;          % dilation radius (px) to close edge gaps before fill
sigmaBlur        = 20;         % Gaussian blur sigma (px)
closeRadius      = 25;         % morphological close disk radius (px)
minArea          = 5000;       % minimum component area (px^2)

% --- Angle x depth grid ---
nAngleBins   = 120;     % angular bins around centroid (covers 0..2*pi)
depthStep_um = 0.5;     % depth bin step (microns)
maxDepth_um  = 20;      % how far inward from cortex to sample (microns)

% --- Selection denoising ---
% Optional Gaussian smoothing along the angle dimension of the per-frame
% (angle x depth) map before picking the max-peak angle. Stabilizes the
% selection against single-bin noise. 0 disables.
angleSmoothBins = 1;

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
                      + double(imread(fullfile(ds{4}(t).folder, ds{4}(t).name)));

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

% Smoothing window for selection
if angleSmoothBins > 0
    smoothWin = 2 * ceil(3 * angleSmoothBins) + 1;
else
    smoothWin = 0;
end

% Per-frame outputs (no 3D cube — only the chosen 1D profile per frame)
maxInward_nm   = nan(nFrames, nDepth);   % rows = frames; cols = depth
avgInward_nm   = nan(nFrames, nDepth);   % angle-averaged comparison
bestAngleIdx   = nan(nFrames, 1);
bestAngle_rad  = nan(nFrames, 1);
peakDepth_um   = nan(nFrames, 1);
peakValue_nm   = nan(nFrames, 1);
centroidXY     = nan(nFrames, 2);
fitRadius_um   = nan(nFrames, 1);

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

    % Angle-averaged inward profile (for comparison kymograph).
    avgInward_nm(fr, :) = mean(slice, 1, 'omitnan');

    % Optional smoothing along the angle axis to stabilize selection.
    if smoothWin > 0
        slice_sel = smoothdata(slice, 1, 'gaussian', smoothWin, 'omitnan');
    else
        slice_sel = slice;
    end

    % Pick the angle whose inward profile peaks highest (scanning inward
    % == scanning along increasing depth).
    perAnglePeak = max(slice_sel, [], 2, 'omitnan');
    [bestPeak, aBest] = max(perAnglePeak, [], 'omitnan');
    if isnan(bestPeak); continue; end

    bestAngleIdx(fr)    = aBest;
    bestAngle_rad(fr)   = angleCenters(aBest);
    maxInward_nm(fr, :) = slice(aBest, :);          % plot raw nm (not the smoothed copy)

    % Depth at which that profile peaks.
    [pkVal, pkIdx]  = max(slice(aBest, :), [], 'omitnan');
    if ~isnan(pkVal)
        peakValue_nm(fr) = pkVal;
        peakDepth_um(fr) = depthAxis_um(pkIdx);
    end

    if mod(fr, max(1, round(nFrames/10))) == 0
        fprintf('  Processed %d / %d frames (%.0f%%)\n', fr, nFrames, 100*fr/nFrames);
    end
end

elapsed = toc;
fprintf('Done! %.1f sec total (%.2f sec/frame)\n\n', elapsed, elapsed/nFrames);

%% ========================== PLOTS =========================================
fprintf('Generating plots...\n');

% --- Plot 1: kymograph of the chosen inward profile per frame ---
fig1 = figure('Position', [100 100 900 500]);
imagesc(depthAxis_um, time_min, maxInward_nm);
set(gca, 'YDir', 'normal');
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Max-peak inward retardance profile vs time', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (nm)';
cb.Label.FontSize = 11;
exportgraphics(fig1, fullfile(outDir, 'inward_kymograph_max_peak.png'), 'Resolution', 200);
savefig(fig1, fullfile(outDir, 'inward_kymograph_max_peak.fig'));
close(fig1);

% --- Plot 2: angle-averaged kymograph (for direct comparison) ---
fig2 = figure('Position', [100 100 900 500]);
imagesc(depthAxis_um, time_min, avgInward_nm);
set(gca, 'YDir', 'normal');
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Angle-averaged inward retardance profile vs time', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (nm)';
cb.Label.FontSize = 11;
exportgraphics(fig2, fullfile(outDir, 'inward_kymograph_avg.png'), 'Resolution', 200);
savefig(fig2, fullfile(outDir, 'inward_kymograph_avg.fig'));
close(fig2);

% --- Plot 3: chosen direction (wave origin angle) over time ---
fig3 = figure('Position', [100 100 800 400]);
plot(time_min, rad2deg(bestAngle_rad), 'k.-', 'LineWidth', 1.2, 'MarkerSize', 10);
xlabel('Time (min)', 'FontSize', 12);
ylabel('Wave origin angle (deg)', 'FontSize', 12);
title('Selected radial direction (max inward peak) over time', 'FontSize', 13);
ylim([0 360]); grid on;
exportgraphics(fig3, fullfile(outDir, 'wave_origin_angle_vs_time.png'), 'Resolution', 200);
savefig(fig3, fullfile(outDir, 'wave_origin_angle_vs_time.fig'));
close(fig3);

% --- Plot 4: peak depth (band trajectory) over time ---
fig4 = figure('Position', [100 100 800 400]);
plot(time_min, peakDepth_um, 'b.-', 'LineWidth', 1.4, 'MarkerSize', 10);
xlabel('Time (min)', 'FontSize', 12);
ylabel('Peak depth from cortex (\mum)', 'FontSize', 12);
title('Wavefront depth vs time (band trajectory)', 'FontSize', 13);
grid on;
exportgraphics(fig4, fullfile(outDir, 'wave_peak_depth_vs_time.png'), 'Resolution', 200);
savefig(fig4, fullfile(outDir, 'wave_peak_depth_vs_time.fig'));
close(fig4);

%% ========================== SAVE DATA =====================================
results = struct();
results.maxInward_nm        = maxInward_nm;
results.avgInward_nm        = avgInward_nm;
results.bestAngleIdx        = bestAngleIdx;
results.bestAngle_rad       = bestAngle_rad;
results.peakDepth_um        = peakDepth_um;
results.peakValue_nm        = peakValue_nm;
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
results.angleSmoothBins     = angleSmoothBins;
results.thresholdMode       = thresholdMode;
results.edgeMethod          = edgeMethod;
results.inputMode           = inputMode;

save(fullfile(outDir, 'wavefront_kymograph_results.mat'), '-struct', 'results');
fprintf('Saved results to %s\n', fullfile(outDir, 'wavefront_kymograph_results.mat'));

% wavefront_kymograph.m
% Per-frame max-peak inward retardance profile -> depth-vs-time kymograph.
%
% Why this is a separate script from contour_retardance.m:
%   contour_retardance.m angle-averages the inward retardance profile, which
%   washes out the localized, weak alignment band that propagates inward
%   during the contraction wave. This script keeps the angle resolution per
%   frame, then uses an angular-anomaly criterion (subtract the per-frame
%   angular mean before selection) so the rotationally symmetric cortex band
%   cancels and the chosen direction is the one carrying the localized
%   wavefront.
%
% The map per frame is computed in one vectorized accumarray pass over the
% interior pixels using the distance transform — no per-angle ray loop, no
% per-boundary-point loop. Only the chosen 1D profile (raw + anomaly) is
% kept per frame (no 3D cube).
%
% Inputs are the same image stacks consumed by contour_retardance.m
% (4-state Polscope is the default; retardance and multipage are also
% supported).
%
% Selection rule: per frame, subtract the angular mean of the (angle x depth)
% map (this kills the rotationally symmetric cortex band), then pick the
% angle whose anomaly profile peaks the highest. The cortex band cancels in
% the anomaly so the chosen direction is driven by the localized wavefront,
% not by the brightest cortex.
%
% Outputs (in outDir):
%   inward_kymograph_chosen_angle_raw.png/.fig      chosen-angle raw inward profile per frame (nm)
%   inward_kymograph_chosen_angle_anomaly.png/.fig  same but with the angular mean subtracted (wavefront only)
%   wave_origin_angle_vs_time.png/.fig              chosen direction over time
%   wave_peak_depth_vs_time.png/.fig                depth where the anomaly peaks vs time -> wave speed
%   wavefront_kymograph_results.mat                 per-frame profiles + selections
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
dt_sec    = 20;
px_per_um = 6.25/2;
retardance_ceiling_nm = 50;
bit_depth = 16;

% --- Segmentation (defaults match contour_retardance.m's working setup) ---
thresholdMode      = 'adaptive';   % 'adaptive', 'edge', 'gradient', 'otsu'
adaptSensitivity   = 0.7;          % adaptthresh sensitivity [0..1]
adaptNeighborhood  = 201;          % adaptthresh neighborhood (odd integer, px)
edgeMethod         = 'Sobel';      % edge() method: 'Canny', 'Sobel', 'Prewitt', 'log'
edgeDilateRadius   = 2;            % dilation radius (px)
gradientPercentile = 70;           % gradient magnitude percentile for 'gradient'
sigmaBlur          = 1;            % Gaussian blur sigma (px)
closeRadius        = 1;            % morphological close disk radius (px)
minArea            = 5000;         % minimum component area (px^2)

% --- Angle x depth grid ---
nAngleBins   = 120;     % angular bins around centroid (covers 0..2*pi)
depthStep_um = 0.5;     % depth bin step (microns)
maxDepth_um  = 10;      % how far inward from cortex to sample (microns)

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

% Smoothing window for selection
if angleSmoothBins > 0
    smoothWin = 2 * ceil(3 * angleSmoothBins) + 1;
else
    smoothWin = 0;
end

% Per-frame outputs (no 3D cube — only the chosen 1D profile per frame)
maxInward_nm     = nan(nFrames, nDepth);   % chosen-angle raw inward profile
anomalyInward_nm = nan(nFrames, nDepth);   % chosen-angle anomaly (slice - per-frame angular mean)
avgInward_nm     = nan(nFrames, nDepth);   % per-frame angular mean (the rotationally symmetric part)
bestAngleIdx     = nan(nFrames, 1);
bestAngle_rad    = nan(nFrames, 1);
peakDepth_um     = nan(nFrames, 1);        % depth where the anomaly peaks (= wavefront depth)
peakValue_nm     = nan(nFrames, 1);        % anomaly value at peakDepth_um
centroidXY       = nan(nFrames, 2);
fitRadius_um     = nan(nFrames, 1);

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

    % --- Segmentation (mirrors contour_retardance.m) ---
    I_blur = imgaussfilt(Iseg, sigmaBlur); %#ok<NASGU>
    I_norm = I_blur / max(I_blur(:));      %#ok<NASGU>

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

        case 'otsu'
            Totsu = graythresh(I_norm);
            if segFromMask
                BW = I_norm < Totsu;
            else
                BW = I_norm > Totsu;
            end

        otherwise
            error('Unknown thresholdMode: %s', thresholdMode);
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

    % Per-frame angular mean (the rotationally symmetric part — cortex band
    % lives mostly here). Subtracting it yields the angular anomaly, which
    % isolates the localized wavefront from the dominant cortex envelope.
    angMean   = mean(slice, 1, 'omitnan');           % 1 x nDepth
    anomaly   = slice - angMean;                     % nA x nDepth
    avgInward_nm(fr, :) = angMean;

    % Optional smoothing along the angle axis to stabilize selection.
    if smoothWin > 0
        anomaly_sel = smoothdata(anomaly, 1, 'gaussian', smoothWin, 'omitnan');
    else
        anomaly_sel = anomaly;
    end

    % Selection: pick the angle whose ANOMALY peaks highest. This rejects
    % the cortex band (which cancels in `anomaly`) so the chosen direction
    % is the one with the strongest localized inward feature.
    perAnglePeak = max(anomaly_sel, [], 2, 'omitnan');
    [bestPeak, aBest] = max(perAnglePeak, [], 'omitnan');
    if isnan(bestPeak); continue; end

    bestAngleIdx(fr)        = aBest;
    bestAngle_rad(fr)       = angleCenters(aBest);
    maxInward_nm(fr, :)     = slice(aBest, :);     % raw nm at the chosen angle
    anomalyInward_nm(fr, :) = anomaly(aBest, :);   % wavefront-only profile at the chosen angle

    % Depth at which the anomaly (wavefront) peaks in the chosen direction.
    [pkVal, pkIdx]  = max(anomaly(aBest, :), [], 'omitnan');
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

% --- Plot 1: raw inward profile at the chosen angle (interpretable nm) ---
% Selection used the angular anomaly so the cortex band doesn't dominate
% the choice; here we plot the actual retardance values for that direction.
fig1 = figure('Position', [100 100 900 500]);
imagesc(depthAxis_um, time_min, maxInward_nm);
set(gca, 'YDir', 'normal');
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Inward retardance profile at chosen angle (raw nm)', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (nm)';
cb.Label.FontSize = 11;
exportgraphics(fig1, fullfile(outDir, 'inward_kymograph_chosen_angle_raw.png'), 'Resolution', 200);
savefig(fig1, fullfile(outDir, 'inward_kymograph_chosen_angle_raw.fig'));
close(fig1);

% --- Plot 2: anomaly inward profile at the chosen angle (wavefront only) ---
% Per-frame angular mean has been subtracted, so the rotationally symmetric
% cortex band is gone and only the localized wavefront remains.
fig2 = figure('Position', [100 100 900 500]);
imagesc(depthAxis_um, time_min, anomalyInward_nm);
set(gca, 'YDir', 'normal');
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Inward retardance anomaly at chosen angle (wavefront only)', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = '\Delta Retardance vs angular mean (nm)';
cb.Label.FontSize = 11;
exportgraphics(fig2, fullfile(outDir, 'inward_kymograph_chosen_angle_anomaly.png'), 'Resolution', 200);
savefig(fig2, fullfile(outDir, 'inward_kymograph_chosen_angle_anomaly.fig'));
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
results.maxInward_nm        = maxInward_nm;          % chosen-angle raw inward profile per frame
results.anomalyInward_nm    = anomalyInward_nm;      % chosen-angle (slice - per-frame angular mean)
results.avgInward_nm        = avgInward_nm;          % per-frame angular mean (the subtracted baseline)
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

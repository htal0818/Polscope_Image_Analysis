% batch_distribution.m
% Batch retardance histogram across many oocytes from static PolScope snapshots.
%
% This script:
%   1. Scans a parent directory for all subfolders containing 'SM' in the name
%   2. Inside each SM folder, looks for a 'Pos0' subfolder
%   3. Loads the first retardance image from Pos0
%   4. Detects the oocyte boundary (same approach as contour_retardance.m)
%   5. Extracts retardance values within the oocyte mask
%   6. Pools all oocyte retardance values and produces a histogram of
%      retardance (nm) vs counts across the entire dataset
%
% Expected folder structure:
%   parent_dir/
%     ├── <something>_SM_<something>/Pos0/*Retardance*
%     ├── <something>_SM_<something>/Pos0/*Retardance*
%     └── ...
%
% REQUIREMENTS:
%   - Image Processing Toolbox
%   - circfit.m (included in this repository)

clear all; close all; clc

%% ========================== USER INPUTS ==================================

% --- Parent directory containing all SM subfolders ---
parent_dir = '/path/to/your/data/';

% --- File pattern for retardance images inside Pos0 ---
retardance_pattern = '*Retardance*';

% --- Retardance calibration ---
retardance_ceiling_nm = 50;   % Polscope retardance ceiling (nm)
bit_depth = 16;               % image bit depth (16-bit = 0..65535)

% --- Spatial calibration ---
px_per_um = 6.25;             % pixels per micron (adjust for your objective)

% --- Boundary detection parameters (from contour_retardance.m) ---
sigmaBlur   = 20;             % Gaussian blur sigma (px) for segmentation
closeRadius = 25;             % morphological close disk radius (px)
minArea     = 5000;           % minimum object area (px^2) to reject debris

% --- Threshold mode ---
% 'otsu'       : automatic Otsu threshold on blurred image (default)
% 'fixed'      : fixed intensity threshold on blurred image
% 'percentile' : threshold at a percentile of blurred image intensity
thresholdMode       = 'otsu';
fixedThreshold      = 500;
percentileThreshold = 30;

% --- Histogram parameters ---
nBins       = 100;            % number of histogram bins
maxRet_nm   = 20;             % max retardance for histogram x-axis (nm)

% --- Output ---
outDir = fullfile(parent_dir, 'batch_retardance_distribution');

% ============================================================================
%% ========================== FIND ALL SM/Pos0 FOLDERS =====================

allDirs = dir(parent_dir);
allDirs = allDirs([allDirs.isdir]);
allDirs = allDirs(~ismember({allDirs.name}, {'.', '..'}));

% Filter to folders containing 'SM' (case-insensitive)
smMask = cellfun(@(n) ~isempty(regexpi(n, 'SM', 'once')), {allDirs.name});
smDirs = allDirs(smMask);

if isempty(smDirs)
    error('No folders containing ''SM'' found in:\n  %s', parent_dir);
end

fprintf('Found %d SM folders in:\n  %s\n\n', numel(smDirs), parent_dir);

%% ========================== SETUP ========================================
if ~exist(outDir, 'dir'); mkdir(outDir); end

maxPixVal = 2^bit_depth - 1;  % 65535 for 16-bit
um_per_px = 1 / px_per_um;

% Pooled retardance values across all oocytes
allRetardance = [];

% Per-oocyte summary statistics
oocyteNames  = {};
oocyteMean   = [];
oocyteMedian = [];
oocyteStd    = [];
oocyteMax    = [];
oocyteArea_um2 = [];

nProcessed = 0;
nSkipped   = 0;

%% ========================== MAIN LOOP ====================================
fprintf('Processing oocytes...\n');
tic;

for si = 1:numel(smDirs)
    smName  = smDirs(si).name;
    pos0Dir = fullfile(parent_dir, smName, 'Pos0');

    % Check that Pos0 exists
    if ~exist(pos0Dir, 'dir')
        fprintf('  [SKIP] %s — no Pos0 folder found\n', smName);
        nSkipped = nSkipped + 1;
        continue;
    end

    % Find retardance image(s) in Pos0
    d = dir(fullfile(pos0Dir, retardance_pattern));
    if isempty(d)
        fprintf('  [SKIP] %s/Pos0 — no retardance images matching "%s"\n', ...
            smName, retardance_pattern);
        nSkipped = nSkipped + 1;
        continue;
    end

    % Sort and take the first retardance image (static snapshot)
    [~, sortIdx] = sort({d.name});
    d = d(sortIdx);
    imgPath = fullfile(d(1).folder, d(1).name);

    %% ----- Load image and convert to retardance (nm) -----
    Iraw = double(imread(imgPath));
    [H, W] = size(Iraw);

    % Convert raw pixel values to retardance in nm
    Iret = (Iraw / maxPixVal) * retardance_ceiling_nm;

    %% ----- Boundary detection (from contour_retardance.m) -----
    switch thresholdMode
        case 'otsu'
            I_blur = imgaussfilt(Iraw, sigmaBlur);
            I_norm = I_blur / max(I_blur(:));
            Totsu  = graythresh(I_norm);
            BW     = I_norm > Totsu;

        case 'fixed'
            I_blur = imgaussfilt(Iraw, sigmaBlur);
            BW     = I_blur > fixedThreshold;

        case 'percentile'
            I_blur = imgaussfilt(Iraw, sigmaBlur);
            pVal   = prctile(I_blur(:), percentileThreshold);
            BW     = I_blur > pVal;

        otherwise
            error('Unknown thresholdMode: %s', thresholdMode);
    end

    % Morphological cleanup
    se = strel('disk', closeRadius);
    BW = imclose(BW, se);
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, minArea);

    % Fallback: gradient-based if Otsu fails
    if ~any(BW(:))
        [Gmag, ~] = imgradient(I_blur);
        thrG = max(2*mean(Gmag(:)), prctile(Gmag(:), 80));
        BW = Gmag >= thrG;
        BW = imclose(BW, se);
        BW = imfill(BW, 'holes');
        BW = bwareaopen(BW, minArea);
    end

    % Keep largest connected component
    L = bwlabel(BW, 8);
    if max(L(:)) >= 1
        S = regionprops(L, 'Area');
        [~, iMax] = max([S.Area]);
        BW = (L == iMax);
    else
        fprintf('  [SKIP] %s — boundary detection failed\n', smName);
        nSkipped = nSkipped + 1;
        continue;
    end

    %% ----- Extract retardance values within the oocyte -----
    retVals = Iret(BW);

    % Pool into the combined dataset
    allRetardance = [allRetardance; retVals(:)];

    % Per-oocyte statistics
    nProcessed = nProcessed + 1;
    oocyteNames{nProcessed}    = smName;
    oocyteMean(nProcessed)     = mean(retVals);
    oocyteMedian(nProcessed)   = median(retVals);
    oocyteStd(nProcessed)      = std(retVals);
    oocyteMax(nProcessed)      = max(retVals);
    oocyteArea_um2(nProcessed) = sum(BW(:)) * um_per_px^2;

    fprintf('  [OK]   %s — %s — mean=%.2f nm, median=%.2f nm, %d px\n', ...
        smName, d(1).name, oocyteMean(nProcessed), oocyteMedian(nProcessed), numel(retVals));
end

elapsed = toc;
fprintf('\nDone! Processed %d oocytes, skipped %d (%.1f sec)\n\n', ...
    nProcessed, nSkipped, elapsed);

if nProcessed == 0
    error('No oocytes were successfully processed. Check your parent_dir and folder structure.');
end

%% ========================== HISTOGRAM: Retardance vs Counts ==============
binEdges = linspace(0, maxRet_nm, nBins + 1);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

fig1 = figure('Position', [100 100 900 550]);
histogram(allRetardance, binEdges, ...
    'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
xlabel('Retardance (nm)', 'FontSize', 13);
ylabel('Counts', 'FontSize', 13);
title(sprintf('Retardance Distribution — %d Oocytes', nProcessed), 'FontSize', 15);
set(gca, 'FontSize', 11);
grid on; box on;

% Add summary stats annotation
grandMean   = mean(allRetardance);
grandMedian = median(allRetardance);
grandStd    = std(allRetardance);
annotation('textbox', [0.62, 0.72, 0.25, 0.15], ...
    'String', sprintf('n = %d oocytes\nMean = %.2f nm\nMedian = %.2f nm\nSD = %.2f nm', ...
        nProcessed, grandMean, grandMedian, grandStd), ...
    'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'FitBoxToText', 'on');

exportgraphics(fig1, fullfile(outDir, 'retardance_histogram_all.png'), 'Resolution', 200);
savefig(fig1, fullfile(outDir, 'retardance_histogram_all.fig'));
close(fig1);

%% ========================== PER-OOCYTE MEAN HISTOGRAM ====================
fig2 = figure('Position', [100 100 900 550]);
histogram(oocyteMean, 20, ...
    'FaceColor', [0.8 0.3 0.2], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
xlabel('Mean Retardance per Oocyte (nm)', 'FontSize', 13);
ylabel('Number of Oocytes', 'FontSize', 13);
title(sprintf('Distribution of Mean Retardance — %d Oocytes', nProcessed), 'FontSize', 15);
set(gca, 'FontSize', 11);
grid on; box on;

exportgraphics(fig2, fullfile(outDir, 'retardance_histogram_per_oocyte_mean.png'), 'Resolution', 200);
savefig(fig2, fullfile(outDir, 'retardance_histogram_per_oocyte_mean.fig'));
close(fig2);

%% ========================== OVERLAY: PER-OOCYTE DISTRIBUTIONS ============
fig3 = figure('Position', [100 100 900 550]);
cmap = parula(nProcessed);
hold on;

for oi = 1:nProcessed
    pos0Dir = fullfile(parent_dir, oocyteNames{oi}, 'Pos0');
    d = dir(fullfile(pos0Dir, retardance_pattern));
    [~, sortIdx] = sort({d.name});
    d = d(sortIdx);
    imgPath = fullfile(d(1).folder, d(1).name);
    Iraw = double(imread(imgPath));
    Iret_i = (Iraw / maxPixVal) * retardance_ceiling_nm;

    % Re-segment to get mask
    switch thresholdMode
        case 'otsu'
            I_blur = imgaussfilt(Iraw, sigmaBlur);
            I_norm = I_blur / max(I_blur(:));
            BW = I_norm > graythresh(I_norm);
        case 'fixed'
            I_blur = imgaussfilt(Iraw, sigmaBlur);
            BW = I_blur > fixedThreshold;
        case 'percentile'
            I_blur = imgaussfilt(Iraw, sigmaBlur);
            BW = I_blur > prctile(I_blur(:), percentileThreshold);
    end
    se = strel('disk', closeRadius);
    BW = imclose(BW, se);
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, minArea);
    L = bwlabel(BW, 8);
    if max(L(:)) >= 1
        S = regionprops(L, 'Area');
        [~, iMax] = max([S.Area]);
        BW = (L == iMax);
    end

    retVals_i = Iret_i(BW);
    [counts_i, ~] = histcounts(retVals_i, binEdges);
    plot(binCenters, counts_i, '-', 'Color', [cmap(oi,:) 0.6], 'LineWidth', 1.2);
end

xlabel('Retardance (nm)', 'FontSize', 13);
ylabel('Counts', 'FontSize', 13);
title(sprintf('Per-Oocyte Retardance Distributions — %d Oocytes', nProcessed), 'FontSize', 15);
set(gca, 'FontSize', 11);
grid on; box on;

exportgraphics(fig3, fullfile(outDir, 'retardance_overlay_per_oocyte.png'), 'Resolution', 200);
savefig(fig3, fullfile(outDir, 'retardance_overlay_per_oocyte.fig'));
close(fig3);

%% ========================== SAVE DATA ====================================
results = struct();
results.allRetardance         = allRetardance;
results.oocyteNames           = oocyteNames;
results.oocyteMean_nm         = oocyteMean;
results.oocyteMedian_nm       = oocyteMedian;
results.oocyteStd_nm          = oocyteStd;
results.oocyteMax_nm          = oocyteMax;
results.oocyteArea_um2        = oocyteArea_um2;
results.nProcessed            = nProcessed;
results.nSkipped              = nSkipped;
results.grandMean_nm          = grandMean;
results.grandMedian_nm        = grandMedian;
results.grandStd_nm           = grandStd;
results.retardance_ceiling_nm = retardance_ceiling_nm;
results.bit_depth             = bit_depth;
results.px_per_um             = px_per_um;
results.binEdges_nm           = binEdges;
results.binCenters_nm         = binCenters;
results.thresholdMode         = thresholdMode;
results.sigmaBlur             = sigmaBlur;
results.closeRadius           = closeRadius;
results.minArea               = minArea;
results.parent_dir            = parent_dir;

save(fullfile(outDir, 'batch_retardance_results.mat'), '-struct', 'results');
fprintf('Saved results to: %s\n', fullfile(outDir, 'batch_retardance_results.mat'));

%% ========================== SUMMARY TABLE ================================
fprintf('\n==================== OOCYTE SUMMARY ====================\n');
fprintf('%-40s %8s %8s %8s %8s\n', 'Oocyte', 'Mean', 'Median', 'SD', 'Max');
fprintf('%-40s %8s %8s %8s %8s\n', '', '(nm)', '(nm)', '(nm)', '(nm)');
fprintf('%s\n', repmat('-', 1, 72));
for oi = 1:nProcessed
    fprintf('%-40s %8.2f %8.2f %8.2f %8.2f\n', ...
        oocyteNames{oi}, oocyteMean(oi), oocyteMedian(oi), ...
        oocyteStd(oi), oocyteMax(oi));
end
fprintf('%s\n', repmat('-', 1, 72));
fprintf('%-40s %8.2f %8.2f %8.2f %8.2f\n', ...
    sprintf('GRAND TOTAL (%d oocytes)', nProcessed), ...
    grandMean, grandMedian, grandStd, max(allRetardance));
fprintf('========================================================\n');
fprintf('\nOutputs saved to: %s\n', outDir);

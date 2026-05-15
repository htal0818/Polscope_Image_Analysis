% batch_distribution.m
% Batch CONTOUR retardance histogram across many oocytes from static PolScope snapshots.
%
% This script:
%   1. Scans a parent directory for all subfolders containing 'SM' in the name
%   2. Inside each SM folder, looks for a 'Pos0' subfolder
%   3. Loads image(s) from Pos0 based on inputMode:
%        'retardance'  — loads a single pre-computed retardance image
%        'four_state'  — loads State1..State4 images, averages them for
%                        segmentation, and loads a separate retardance image
%                        for measurement
%   4. Calls measure_contour_retardance() to detect the oocyte boundary
%      and sample retardance (nm) along the cortical contour
%   5. Pools contour retardance values across all oocytes and produces
%      a histogram of retardance (nm) vs counts
%
% Uses the same segmentation + contour measurement pipeline as
% contour_retardance.m via the shared function measure_contour_retardance().
%
% Expected folder structure:
%   parent_dir/
%     ├── <something>_SM_<something>/Pos0/*Retardance*
%     ├── <something>_SM_<something>/Pos0/*Retardance*
%     └── ...
%
% REQUIREMENTS:
%   - Image Processing Toolbox
%   - measure_contour_retardance.m (in this repository)
%   - circfit.m (in this repository)

clear all; close all; clc

%% ========================== USER INPUTS ==================================

% --- Parent directory containing all SM subfolders ---
parent_dir = '/path/to/your/data/';

% --- Input mode ---
% 'retardance'  : pre-computed retardance images (default, current behavior)
% 'four_state'  : raw 4-state Polscope data — State1..State4 summed/averaged
%                 for segmentation, separate retardance image for measurement
inputMode = 'retardance';

% --- File pattern for retardance images inside Pos0 ---
retardance_pattern = '*Retardance*';

% --- File patterns for four_state mode ---
state_patterns = {'*State1*', '*State2*', '*State3*', '*State4*'};

% --- Options passed to measure_contour_retardance() ---
% (see measure_contour_retardance.m for full list and defaults)
opts = struct();
opts.retardance_ceiling_nm = 50;    % Polscope retardance ceiling (nm)
opts.bit_depth             = 16;    % image bit depth (16-bit = 0..65535)
opts.sigmaBlur             = 20;    % Gaussian blur sigma (px) for segmentation
opts.closeRadius           = 25;    % morphological close disk radius (px)
opts.minArea               = 5000;  % minimum object area (px^2) to reject debris
opts.thresholdMode         = 'otsu';
opts.fixedThreshold        = 500;
opts.percentileThreshold   = 30;
opts.adaptiveSensitivity   = 0.5;  % for 'adaptive' mode (0-1, higher = more foreground)

% --- Spatial calibration ---
px_per_um = 6.25;             % pixels per micron (adjust for your objective)

% --- Cortical inset depth ---
boundaryInset_um           = 2.5;   % cortical depth (um), typically 2-3 um
opts.boundaryInset_px      = round(boundaryInset_um * px_per_um);

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

um_per_px = 1 / px_per_um;

% Pooled contour retardance values across all oocytes
allContourRet = [];

% Per-oocyte summary statistics
oocyteNames    = {};
oocyteMean     = [];
oocyteMedian   = [];
oocyteStd      = [];
oocyteMax      = [];
oocyteMin      = [];
oocyteR_um     = [];   % circle-fit radius per oocyte

% Store per-oocyte contour values for overlay plot
perOocyteContour = {};

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

    % Load image(s) based on input mode
    imgLabel = '';
    switch inputMode
        case 'retardance'
            % Find retardance image(s) in Pos0
            d = dir(fullfile(pos0Dir, retardance_pattern));
            d = d(~[d.isdir]);  % exclude directories
            if isempty(d)
                fprintf('  [SKIP] %s/Pos0 — no retardance images matching "%s"\n', ...
                    smName, retardance_pattern);
                nSkipped = nSkipped + 1;
                continue;
            end
            [~, sortIdx] = sort({d.name});
            d = d(sortIdx);
            imgPath = fullfile(d(1).folder, d(1).name);
            Iraw = double(imread(imgPath));
            imgLabel = d(1).name;

        case 'four_state'
            % Find State1..State4 images, sum and average for segmentation
            stateFiles = cell(1, 4);
            stateMissing = false;
            for qi = 1:4
                ds = dir(fullfile(pos0Dir, state_patterns{qi}));
                ds = ds(~[ds.isdir]);  % exclude directories
                if isempty(ds)
                    fprintf('  [SKIP] %s/Pos0 — no images matching "%s"\n', ...
                        smName, state_patterns{qi});
                    stateMissing = true;
                    break;
                end
                [~, idx] = sort({ds.name});
                ds = ds(idx);
                stateFiles{qi} = fullfile(ds(1).folder, ds(1).name);
            end
            if stateMissing
                nSkipped = nSkipped + 1;
                continue;
            end
            Isum = (double(imread(stateFiles{1})) ...
                  + double(imread(stateFiles{2})) ...
                  + double(imread(stateFiles{3})) ...
                  + double(imread(stateFiles{4}))) / 4;
            opts.Iseg = Isum;
            opts.thresholdMode = 'adaptive';

            % Also load the retardance image for measurement
            d = dir(fullfile(pos0Dir, retardance_pattern));
            d = d(~[d.isdir]);  % exclude directories
            if isempty(d)
                fprintf('  [SKIP] %s/Pos0 — no retardance images matching "%s"\n', ...
                    smName, retardance_pattern);
                nSkipped = nSkipped + 1;
                continue;
            end
            [~, sortIdx] = sort({d.name});
            d = d(sortIdx);
            Iraw = double(imread(fullfile(d(1).folder, d(1).name)));
            imgLabel = sprintf('%s (four_state seg)', d(1).name);

        otherwise
            error('Unknown inputMode: %s. Use ''retardance'' or ''four_state''.', inputMode);
    end

    % Call the shared contour measurement function
    res = measure_contour_retardance(Iraw, opts);

    if ~res.success
        fprintf('  [SKIP] %s — boundary detection failed\n', smName);
        nSkipped = nSkipped + 1;
        continue;
    end

    % Pool contour values into the combined dataset
    cv = res.contourValues;
    allContourRet = [allContourRet; cv(:)];

    % Per-oocyte statistics
    nProcessed = nProcessed + 1;
    oocyteNames{nProcessed}      = smName;
    oocyteMean(nProcessed)        = res.contourMean;
    oocyteMedian(nProcessed)      = median(cv);
    oocyteStd(nProcessed)         = res.contourStd;
    oocyteMax(nProcessed)         = res.contourMax;
    oocyteMin(nProcessed)         = res.contourMin;
    oocyteR_um(nProcessed)        = res.R_fit * um_per_px;
    perOocyteContour{nProcessed}  = cv;

    fprintf('  [OK]   %s — %s — contour mean=%.2f nm, R=%.0f um, %d pts\n', ...
        smName, imgLabel, res.contourMean, oocyteR_um(nProcessed), numel(cv));
end

elapsed = toc;
fprintf('\nDone! Processed %d oocytes, skipped %d (%.1f sec)\n\n', ...
    nProcessed, nSkipped, elapsed);

if nProcessed == 0
    error('No oocytes were successfully processed. Check your parent_dir and folder structure.');
end

%% ========================== HISTOGRAM: Contour Retardance vs Counts ======
binEdges   = linspace(0, maxRet_nm, nBins + 1);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

grandMean   = mean(allContourRet);
grandMedian = median(allContourRet);
grandStd    = std(allContourRet);

fig1 = figure('Position', [100 100 900 550]);
histogram(allContourRet, binEdges, ...
    'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
xlabel('Contour Retardance (nm)', 'FontSize', 13);
ylabel('Counts', 'FontSize', 13);
title(sprintf('Contour Retardance Distribution — %d Oocytes', nProcessed), 'FontSize', 15);
set(gca, 'FontSize', 11);
grid on; box on;

annotation('textbox', [0.62, 0.72, 0.25, 0.15], ...
    'String', sprintf('n = %d oocytes\nMean = %.2f nm\nMedian = %.2f nm\nSD = %.2f nm', ...
        nProcessed, grandMean, grandMedian, grandStd), ...
    'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'FitBoxToText', 'on');

exportgraphics(fig1, fullfile(outDir, 'contour_retardance_histogram_all.png'), 'Resolution', 200);
savefig(fig1, fullfile(outDir, 'contour_retardance_histogram_all.fig'));
close(fig1);

%% ========================== PER-OOCYTE MEAN HISTOGRAM ====================
fig2 = figure('Position', [100 100 900 550]);
histogram(oocyteMean, 20, ...
    'FaceColor', [0.8 0.3 0.2], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
xlabel('Mean Contour Retardance per Oocyte (nm)', 'FontSize', 13);
ylabel('Number of Oocytes', 'FontSize', 13);
title(sprintf('Distribution of Mean Contour Retardance — %d Oocytes', nProcessed), 'FontSize', 15);
set(gca, 'FontSize', 11);
grid on; box on;

exportgraphics(fig2, fullfile(outDir, 'contour_retardance_histogram_per_oocyte_mean.png'), 'Resolution', 200);
savefig(fig2, fullfile(outDir, 'contour_retardance_histogram_per_oocyte_mean.fig'));
close(fig2);

%% ========================== PLOT 3: MEAN + INDIVIDUAL TRACES =============
% Each oocyte as a thin gray line, population mean as a bold colored line
fig3 = figure('Position', [100 100 900 550]);
hold on;

% Normalize each oocyte to probability density so all eggs contribute equally
perOocyteDensity = zeros(nProcessed, nBins);
for oi = 1:nProcessed
    [counts_i, ~] = histcounts(perOocyteContour{oi}, binEdges, 'Normalization', 'probability');
    perOocyteDensity(oi,:) = counts_i;
    plot(binCenters, counts_i, '-', 'Color', [0.6 0.6 0.6 0.35], 'LineWidth', 0.8);
end

% Population mean + SD band
meanDensity = mean(perOocyteDensity, 1);
sdDensity   = std(perOocyteDensity, 0, 1);

fill([binCenters fliplr(binCenters)], ...
     [meanDensity + sdDensity, fliplr(meanDensity - sdDensity)], ...
     [0.2 0.4 0.8], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(binCenters, meanDensity, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);

xlabel('Contour Retardance (nm)', 'FontSize', 13);
ylabel('Fraction of Contour', 'FontSize', 13);
title(sprintf('Per-Oocyte Contour Distributions — %d Oocytes', nProcessed), 'FontSize', 15);
legend({'Individual oocytes', '\pm1 SD', 'Population mean'}, ...
    'Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11);
grid on; box on;

exportgraphics(fig3, fullfile(outDir, 'contour_retardance_mean_traces.png'), 'Resolution', 200);
savefig(fig3, fullfile(outDir, 'contour_retardance_mean_traces.fig'));
close(fig3);

%% ========================== PLOT 4: HEATMAP (oocyte x retardance) ========
% Each row is one oocyte's contour distribution, sorted by mean retardance
fig4 = figure('Position', [100 100 900 600]);

[~, sortOrder] = sort(oocyteMean);
sortedDensity  = perOocyteDensity(sortOrder, :);
sortedNames    = oocyteNames(sortOrder);

imagesc(binCenters, 1:nProcessed, sortedDensity);
set(gca, 'YDir', 'normal');
colormap parula; cb = colorbar;
cb.Label.String = 'Fraction of Contour';
cb.Label.FontSize = 11;

xlabel('Contour Retardance (nm)', 'FontSize', 13);
ylabel('Oocyte (sorted by mean retardance)', 'FontSize', 13);
title(sprintf('Contour Retardance Heatmap — %d Oocytes', nProcessed), 'FontSize', 15);
set(gca, 'FontSize', 11);

% Label y-axis with oocyte names if <= 30, otherwise use numbers
if nProcessed <= 30
    set(gca, 'YTick', 1:nProcessed, 'YTickLabel', sortedNames, 'FontSize', 7);
end

exportgraphics(fig4, fullfile(outDir, 'contour_retardance_heatmap.png'), 'Resolution', 200);
savefig(fig4, fullfile(outDir, 'contour_retardance_heatmap.fig'));
close(fig4);

%% ========================== SAVE DATA ====================================
results = struct();
results.allContourRet         = allContourRet;
results.perOocyteContour      = {perOocyteContour};
results.oocyteNames           = oocyteNames;
results.oocyteMean_nm         = oocyteMean;
results.oocyteMedian_nm       = oocyteMedian;
results.oocyteStd_nm          = oocyteStd;
results.oocyteMax_nm          = oocyteMax;
results.oocyteMin_nm          = oocyteMin;
results.oocyteR_um            = oocyteR_um;
results.nProcessed            = nProcessed;
results.nSkipped              = nSkipped;
results.grandMean_nm          = grandMean;
results.grandMedian_nm        = grandMedian;
results.grandStd_nm           = grandStd;
results.retardance_ceiling_nm = opts.retardance_ceiling_nm;
results.bit_depth             = opts.bit_depth;
results.px_per_um             = px_per_um;
results.binEdges_nm           = binEdges;
results.binCenters_nm         = binCenters;
results.opts                  = opts;
results.parent_dir            = parent_dir;

save(fullfile(outDir, 'batch_retardance_results.mat'), '-struct', 'results');
fprintf('Saved results to: %s\n', fullfile(outDir, 'batch_retardance_results.mat'));

%% ========================== SUMMARY TABLE ================================
fprintf('\n==================== OOCYTE SUMMARY ====================\n');
fprintf('%-40s %8s %8s %8s %8s %8s\n', 'Oocyte', 'Mean', 'Median', 'SD', 'Max', 'R(um)');
fprintf('%-40s %8s %8s %8s %8s %8s\n', '', '(nm)', '(nm)', '(nm)', '(nm)', '');
fprintf('%s\n', repmat('-', 1, 80));
for oi = 1:nProcessed
    fprintf('%-40s %8.2f %8.2f %8.2f %8.2f %8.1f\n', ...
        oocyteNames{oi}, oocyteMean(oi), oocyteMedian(oi), ...
        oocyteStd(oi), oocyteMax(oi), oocyteR_um(oi));
end
fprintf('%s\n', repmat('-', 1, 80));
fprintf('%-40s %8.2f %8.2f %8.2f %8.2f\n', ...
    sprintf('GRAND TOTAL (%d oocytes)', nProcessed), ...
    grandMean, grandMedian, grandStd, max(allContourRet));
fprintf('========================================================\n');
fprintf('\nOutputs saved to: %s\n', outDir);

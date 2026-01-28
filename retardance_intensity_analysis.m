% retardance_intensity_analysis.m
% Analyze raw intensity statistics from Polscope retardance images over time.
%
% This script:
% 1. Reads Polscope image stacks with '1_Retardance' in filename
% 2. Calculates mean intensity and standard deviation for each frame
% 3. Plots mean intensity and std deviation over time
%
% Author: Generated for Polscope Image Analysis pipeline
% Date: January 2026

clear all; close all; clc

%% ========================== USER INPUTS ==================================
% --- Directory containing Polscope images ---
base_dir = '/Users/hridaytalreja/Desktop/Jan_data_2026/jan_20_2026_FSW_and_eggs_50msexp_15sint_20x_50nmceiling/eggs/SMS_2026_0120_1518_1/Pos0/';

% --- Search pattern for retardance images ---
retardance_pattern = '*1_Retardance*';

% --- Timing parameters ---
dt_sec = 15;       % seconds per frame (adjust to your acquisition rate)

% --- Optional: crop region of interest (set doCrop=false to use full image) ---
doCrop = false;
cropRect = [500 500 2000 2000];  % [x y w h] in pixels

% --- Output directory ---
outDir = fullfile(base_dir, 'retardance_intensity_analysis');

% ============================================================================
%% Find retardance images
fprintf('Searching for retardance images in:\n  %s\n', base_dir);
search_pattern = fullfile(base_dir, retardance_pattern);
d_ret = dir(search_pattern);

nFrames = numel(d_ret);
fprintf('Found %d retardance images matching pattern "%s"\n\n', nFrames, retardance_pattern);

if nFrames == 0
    error('No images found! Check base_dir and retardance_pattern.');
end

% Sort files by name to ensure correct temporal order
[~, sortIdx] = sort({d_ret.name});
d_ret = d_ret(sortIdx);

%% Preallocate arrays for statistics
mean_intensity = nan(nFrames, 1);
std_intensity  = nan(nFrames, 1);
time_sec       = (0:nFrames-1)' * dt_sec;
time_min       = time_sec / 60;

%% Create output directory
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('Created output directory: %s\n', outDir);
end

%% Process each frame
fprintf('Processing %d frames...\n', nFrames);
tic;

for fr = 1:nFrames
    % Load image
    imgPath = fullfile(d_ret(fr).folder, d_ret(fr).name);
    img = double(imread(imgPath));

    % Optional crop
    if doCrop
        img = imcrop(img, cropRect);
    end

    % Calculate statistics on raw intensity
    mean_intensity(fr) = mean(img(:));
    std_intensity(fr)  = std(img(:));

    % Progress update
    if mod(fr, 50) == 0 || fr == nFrames
        fprintf('  Processed %d / %d frames (%.1f%%)\n', fr, nFrames, 100*fr/nFrames);
    end
end

elapsed = toc;
fprintf('Done! Processing time: %.1f seconds (%.2f sec/frame)\n\n', elapsed, elapsed/nFrames);

%% Create plots
fprintf('Generating plots...\n');

% Figure 1: Mean Intensity over Time
fig1 = figure('Name', 'Mean Retardance Intensity', 'Position', [100 100 800 500]);
plot(time_min, mean_intensity, 'b-', 'LineWidth', 1.5);
hold on;
plot(time_min, mean_intensity, 'b.', 'MarkerSize', 8);
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Mean Intensity (a.u.)', 'FontSize', 12);
title('Mean Retardance Intensity Over Time', 'FontSize', 14);
grid on;
xlim([0 max(time_min)]);

% Add statistics annotation
stats_text = sprintf('Mean: %.2f\nMin: %.2f\nMax: %.2f', ...
    mean(mean_intensity), min(mean_intensity), max(mean_intensity));
text(0.02, 0.98, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

% Save figure 1
saveas(fig1, fullfile(outDir, 'mean_intensity_over_time.png'));
saveas(fig1, fullfile(outDir, 'mean_intensity_over_time.fig'));

% Figure 2: Standard Deviation over Time
fig2 = figure('Name', 'Retardance Intensity Std Dev', 'Position', [150 150 800 500]);
plot(time_min, std_intensity, 'r-', 'LineWidth', 1.5);
hold on;
plot(time_min, std_intensity, 'r.', 'MarkerSize', 8);
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Standard Deviation (a.u.)', 'FontSize', 12);
title('Retardance Intensity Standard Deviation Over Time', 'FontSize', 14);
grid on;
xlim([0 max(time_min)]);

% Add statistics annotation
stats_text = sprintf('Mean Std: %.2f\nMin Std: %.2f\nMax Std: %.2f', ...
    mean(std_intensity), min(std_intensity), max(std_intensity));
text(0.02, 0.98, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

% Save figure 2
saveas(fig2, fullfile(outDir, 'std_intensity_over_time.png'));
saveas(fig2, fullfile(outDir, 'std_intensity_over_time.fig'));

% Figure 3: Combined plot with dual y-axis
fig3 = figure('Name', 'Combined Retardance Statistics', 'Position', [200 200 900 500]);
yyaxis left
plot(time_min, mean_intensity, 'b-', 'LineWidth', 1.5);
ylabel('Mean Intensity (a.u.)', 'FontSize', 12);
yyaxis right
plot(time_min, std_intensity, 'r-', 'LineWidth', 1.5);
ylabel('Standard Deviation (a.u.)', 'FontSize', 12);
xlabel('Time (minutes)', 'FontSize', 12);
title('Retardance Intensity Statistics Over Time', 'FontSize', 14);
legend('Mean Intensity', 'Std Deviation', 'Location', 'best');
grid on;
xlim([0 max(time_min)]);

% Save figure 3
saveas(fig3, fullfile(outDir, 'combined_intensity_stats.png'));
saveas(fig3, fullfile(outDir, 'combined_intensity_stats.fig'));

%% Save data to .mat file
output_data = struct();
output_data.mean_intensity = mean_intensity;
output_data.std_intensity = std_intensity;
output_data.time_sec = time_sec;
output_data.time_min = time_min;
output_data.nFrames = nFrames;
output_data.dt_sec = dt_sec;
output_data.base_dir = base_dir;
output_data.retardance_pattern = retardance_pattern;
output_data.filenames = {d_ret.name};
output_data.doCrop = doCrop;
if doCrop
    output_data.cropRect = cropRect;
end

save(fullfile(outDir, 'retardance_intensity_data.mat'), '-struct', 'output_data');
fprintf('Saved data to: %s\n', fullfile(outDir, 'retardance_intensity_data.mat'));

%% Print summary
fprintf('\n========== SUMMARY ==========\n');
fprintf('Dataset: %s\n', base_dir);
fprintf('Pattern: %s\n', retardance_pattern);
fprintf('Frames processed: %d\n', nFrames);
fprintf('Total duration: %.2f minutes\n', max(time_min));
fprintf('Frame interval: %.1f seconds\n', dt_sec);
fprintf('\nIntensity Statistics:\n');
fprintf('  Mean intensity: %.2f +/- %.2f (range: %.2f - %.2f)\n', ...
    mean(mean_intensity), std(mean_intensity), min(mean_intensity), max(mean_intensity));
fprintf('  Std deviation:  %.2f +/- %.2f (range: %.2f - %.2f)\n', ...
    mean(std_intensity), std(std_intensity), min(std_intensity), max(std_intensity));
fprintf('\nOutputs saved to: %s\n', outDir);
fprintf('=============================\n');

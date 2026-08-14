% nematic_order_analysis.m
% Nematic order parameter analysis from LC-PolScope image stacks.
%
% Computes the 2D nematic Q tensor and scalar order parameter S following
% Mirza et al. (eLife 2024, arXiv:2306.15352) from slow-axis orientation
% and retardance images.
%
% The per-pixel nematic Q tensor is:
%     Q_ij = S_pixel * (n_i n_j - delta_ij / 2)
% where S_pixel = retardance at each pixel.  Independent components:
%     q1 = retardance * cos(2*phi) / 2    (Q_11 = -Q_22)
%     q2 = retardance * sin(2*phi) / 2    (Q_12 =  Q_21)
% The scalar order parameter is:
%     S = sqrt(2 * Q_ij * Q_ij) = 2 * sqrt(q1^2 + q2^2)
%
% For each frame, computes:
%   - Local S map (Gaussian-averaged Q tensor)
%   - Global S (whole-mask average)
%   - Tangential/radial alignment angle
%   - Cortex angular profile S(theta)
%   - Q-tensor components q1, q2
%
% Outputs:
%   - nematic_results.mat with all time-series and spatial data
%   - PNG overlays for example frames
%
% REQUIREMENTS:
%   - Image Processing Toolbox
%   - circfit.m (included in this repository)

clear all; close all; clc

%% ========================== USER INPUTS ==================================

base_dir = '/path/to/your/data/Pos0/';

% --- File discovery patterns ---
sa_pattern  = '*Slow Axis Orientation*';
ret_pattern = '*Retardance*';

% --- Segmentation ---
% 'states'  : use State1-4 images (summed, inverted Otsu)
% 'sa'      : segment from slow-axis image itself
% 'ret'     : segment from retardance image
% 'mask'    : load pre-computed mask from file
segMode = 'states';
mask_file = '';  % path to pre-computed mask (only used if segMode='mask')

% --- Timing & calibration ---
dt_sec    = 15;        % seconds per frame
px_per_um = 3.125;     % pixels per micron

% --- Nematic analysis parameters ---
sigma_um         = 5.0;   % Gaussian sigma for local S map (microns)
n_theta_bins     = 100;   % angular bins around cortex
director_spacing = 40;    % grid spacing for director overlay (pixels)
inset_px         = 10;    % inward shift of boundary sampling (pixels)

% --- Segmentation parameters ---
sigmaBlur   = 20;      % Gaussian blur sigma (px)
closeRadius = 25;      % morphological close disk radius (px)
minArea     = 5000;    % minimum object area (px^2)

% --- Mask caching ---
cacheForceRecalcEveryN = 25;  % re-segment every N frames

% --- Encoding ---
% 'auto', 'openpolscope', 'degrees', 'radians'
encoding = 'auto';

%% ========================== FILE DISCOVERY ================================
fprintf('Base directory: %s\n', base_dir);

sa_files  = dir(fullfile(base_dir, sa_pattern));
sa_files  = sa_files(~[sa_files.isdir]);
% Filter to image files only
imgExt = {'.tif','.tiff','.png','.jpg'};
keep = false(size(sa_files));
for i = 1:numel(sa_files)
    [~,~,ext] = fileparts(sa_files(i).name);
    keep(i) = any(strcmpi(ext, imgExt));
end
sa_files = sa_files(keep);
[~, idx] = sort({sa_files.name});
sa_files = sa_files(idx);

ret_files = dir(fullfile(base_dir, ret_pattern));
ret_files = ret_files(~[ret_files.isdir]);
keep = false(size(ret_files));
for i = 1:numel(ret_files)
    [~,~,ext] = fileparts(ret_files(i).name);
    keep(i) = any(strcmpi(ext, imgExt));
end
ret_files = ret_files(keep);
[~, idx] = sort({ret_files.name});
ret_files = ret_files(idx);

nFrames = numel(sa_files);
fprintf('  Slow-axis files:  %d\n', nFrames);
fprintf('  Retardance files: %d\n', numel(ret_files));

if nFrames == 0
    error('No slow-axis files found matching pattern "%s" in %s', sa_pattern, base_dir);
end

% State images for segmentation
state_patterns = {'*State1*', '*State2*', '*State3*', '*State4*'};
state_files = cell(1,4);
has_states = true;
for s = 1:4
    tmp = dir(fullfile(base_dir, state_patterns{s}));
    tmp = tmp(~[tmp.isdir]);
    if isempty(tmp)
        has_states = false;
    else
        state_files{s} = tmp;
    end
end

%% ========================== DETECT ENCODING ===============================
sa_img0 = imread(fullfile(sa_files(1).folder, sa_files(1).name));
if size(sa_img0,3) > 1, sa_img0 = sa_img0(:,:,1); end
sa_img0 = double(sa_img0);
[H, W] = size(sa_img0);

if strcmp(encoding, 'auto')
    mx = max(sa_img0(:));
    if isfloat(sa_img0) && mx <= 2*pi + 0.1
        scale = 180/pi; enc_name = 'float_radians';
    elseif mx <= 180.5
        scale = 1.0; enc_name = 'float_degrees_180';
    elseif mx > 17500 && mx <= 18001
        scale = 0.01; enc_name = 'uint16_openpolscope';
    elseif mx > 30000
        scale = 180/65535; enc_name = 'uint16_full_range';
    else
        scale = 0.01; enc_name = 'default_openpolscope';
    end
elseif strcmp(encoding, 'openpolscope')
    scale = 0.01; enc_name = 'openpolscope';
elseif strcmp(encoding, 'degrees')
    scale = 1.0; enc_name = 'degrees';
elseif strcmp(encoding, 'radians')
    scale = 180/pi; enc_name = 'radians';
end
fprintf('Frames: %d, size: %dx%d, encoding: %s\n', nFrames, H, W, enc_name);

%% ========================== SETUP =========================================

sigma_px = sigma_um * px_per_um;

out_dir = fullfile(base_dir, 'nematic_analysis_out');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

time_sec = (0:nFrames-1)' * dt_sec;
time_min = time_sec / 60;
theta_centers_deg = linspace(0, 360, n_theta_bins+1);
theta_centers_deg = (theta_centers_deg(1:end-1) + theta_centers_deg(2:end)) / 2;

% Preallocate time-series
S_global_ts    = nan(nFrames, 1);
psi_global_ts  = nan(nFrames, 1);
q1_global_ts   = nan(nFrames, 1);
q2_global_ts   = nan(nFrames, 1);
align_mean_ts  = nan(nFrames, 1);
S_kymo         = nan(nFrames, n_theta_bins);
psi_kymo       = nan(nFrames, n_theta_bins);
q1_kymo        = nan(nFrames, n_theta_bins);
q2_kymo        = nan(nFrames, n_theta_bins);

% Example frames for spatial overlays
example_frames = unique([1, round(nFrames/2), nFrames]);

mask = [];

%% ========================== MAIN LOOP =====================================
fprintf('\nProcessing %d frames...\n', nFrames);
tic;

for fr = 1:nFrames

    % --- Load slow-axis ---
    sa_img = imread(fullfile(sa_files(fr).folder, sa_files(fr).name));
    if size(sa_img,3) > 1, sa_img = sa_img(:,:,1); end
    sa_float = double(sa_img);

    phi_deg = mod(sa_float * scale, 180);
    phi_rad = deg2rad(phi_deg);

    % --- Load retardance ---
    if fr <= numel(ret_files)
        retardance = imread(fullfile(ret_files(fr).folder, ret_files(fr).name));
        if size(retardance,3) > 1, retardance = mean(double(retardance),3); end
        retardance = double(retardance);
        if ~isequal(size(retardance), [H W])
            retardance = imresize(retardance, [H W]);
        end
        weight = retardance;
    else
        warning('Frame %d: no retardance image — falling back to uniform weights.', fr);
        weight = ones(H, W);
    end

    % --- Segmentation ---
    if isempty(mask) || mod(fr-1, cacheForceRecalcEveryN) == 0
        switch segMode
            case 'states'
                if has_states
                    mask = segment_from_states_m(base_dir, state_patterns, ...
                                                 sigmaBlur, closeRadius, minArea);
                else
                    mask = segment_from_image_m(weight, sigmaBlur, closeRadius, minArea);
                end
            case 'sa'
                mask = segment_from_image_m(phi_deg, sigmaBlur, closeRadius, minArea);
            case 'ret'
                mask = segment_from_image_m(weight, sigmaBlur, closeRadius, minArea);
            case 'mask'
                if fr == 1
                    mask = imread(mask_file) > 0;
                    if size(mask,3) > 1, mask = mask(:,:,1) > 0; end
                end
        end
        if ~isequal(size(mask), [H W])
            mask = imresize(double(mask), [H W]) > 0.5;
        end
    end

    if ~any(mask(:))
        if mod(fr-1, 25) == 0
            fprintf('  Frame %d: no mask, skipping\n', fr);
        end
        continue
    end

    % --- Boundary ---
    [xb, yb, xc, yc, R_fit] = find_boundary_and_center_m(mask);
    if isempty(xb), continue; end

    % --- Nematic order (local map) ---
    [S_map, psi_map, q1_map, q2_map] = nematic_order_map_m(phi_rad, weight, mask, sigma_px);

    % --- Nematic order (global) ---
    [S_g, psi_g, q1_g, q2_g] = nematic_order_global_m(phi_rad, weight, mask);

    % --- Tangential/radial alignment ---
    alpha_map = tangential_radial_alignment_m(phi_rad, mask, xc, yc);

    % --- Cortex nematic profile ---
    [S_theta, psi_theta, q1_theta, q2_theta] = ...
        cortex_nematic_profile_m(phi_rad, weight, xb, yb, xc, yc, n_theta_bins, inset_px);

    % --- Store time-series ---
    S_global_ts(fr)    = S_g;
    psi_global_ts(fr)  = psi_g;
    q1_global_ts(fr)   = q1_g;
    q2_global_ts(fr)   = q2_g;

    alpha_valid = alpha_map(mask & ~isnan(alpha_map));
    if ~isempty(alpha_valid)
        align_mean_ts(fr) = rad2deg(nanmean(alpha_valid));
    end

    S_kymo(fr,:)   = S_theta;
    psi_kymo(fr,:) = psi_theta;
    q1_kymo(fr,:)  = q1_theta;
    q2_kymo(fr,:)  = q2_theta;

    % --- Save spatial overlays for example frames ---
    if ismember(fr, example_frames)
        frame_dir = fullfile(out_dir, sprintf('frame_%04d', fr-1));
        if ~exist(frame_dir,'dir'), mkdir(frame_dir); end

        save_frame_plots(phi_deg, phi_rad, weight, mask, sa_img, ...
                         xb, yb, xc, yc, R_fit, ...
                         S_map, psi_map, S_g, psi_g, ...
                         alpha_map, S_theta, psi_theta, theta_centers_deg, ...
                         q1_map, q2_map, q1_g, q2_g, q1_theta, q2_theta, ...
                         frame_dir, px_per_um, sigma_um, director_spacing);
    end

    % --- Progress ---
    if mod(fr-1, 25) == 0 || fr == nFrames
        elapsed = toc;
        fprintf('  Frame %d/%d  S_global=%.4f  (%.1fs elapsed)\n', ...
                fr, nFrames, S_g, elapsed);
    end
end

%% ========================== TIME-SERIES PLOTS =============================

if nFrames > 1
    % S_global time series
    figure('Visible','off');
    subplot(2,1,1);
    plot(time_min, S_global_ts, 'b-', 'LineWidth', 1.5);
    xlabel('Time (min)'); ylabel('S_{global}');
    title('Global nematic order parameter');

    subplot(2,1,2);
    plot(time_min, align_mean_ts, 'r-', 'LineWidth', 1.5);
    xlabel('Time (min)'); ylabel('Alignment angle (deg)');
    title('Mean alignment angle (0=radial, 90=tangential)');
    saveas(gcf, fullfile(out_dir, 'timeseries_S_global.png'));
    close;

    % S kymograph
    figure('Visible','off');
    imagesc(theta_centers_deg, time_min, S_kymo);
    colormap hot; colorbar;
    xlabel('Angle (deg)'); ylabel('Time (min)');
    title('S(\theta) kymograph');
    saveas(gcf, fullfile(out_dir, 'kymograph_S.png'));
    close;

    % Psi kymograph
    figure('Visible','off');
    imagesc(theta_centers_deg, time_min, rad2deg(psi_kymo));
    colormap hsv; colorbar;
    xlabel('Angle (deg)'); ylabel('Time (min)');
    title('\psi(\theta) kymograph');
    saveas(gcf, fullfile(out_dir, 'kymograph_psi.png'));
    close;

    fprintf('  Saved: timeseries_S_global.png, kymograph_S.png, kymograph_psi.png\n');
end

%% ========================== SAVE DATA =====================================

save(fullfile(out_dir, 'nematic_results.mat'), ...
     'S_global_ts', 'psi_global_ts', 'q1_global_ts', 'q2_global_ts', ...
     'align_mean_ts', ...
     'S_kymo', 'psi_kymo', 'q1_kymo', 'q2_kymo', ...
     'theta_centers_deg', 'time_sec', 'time_min', ...
     'nFrames', 'dt_sec', 'sigma_um', 'px_per_um', 'enc_name');

%% ========================== SUMMARY =======================================
valid = ~isnan(S_global_ts);
fprintf('\n========== NEMATIC ANALYSIS SUMMARY ==========\n');
fprintf('  Frames processed:      %d / %d\n', sum(valid), nFrames);
if nFrames > 1
    fprintf('  Duration:              %.1f min (dt=%ds)\n', time_min(end), dt_sec);
end
fprintf('  Mean S_global:         %.4f +/- %.4f\n', nanmean(S_global_ts), nanstd(S_global_ts));
fprintf('  Mean alignment angle:  %.1f deg (0=radial, 90=tangential)\n', nanmean(align_mean_ts));
fprintf('  Output:                %s\n', out_dir);
fprintf('================================================\n');
fprintf('  Saved: nematic_results.mat\n');


%% ========================== LOCAL FUNCTIONS ================================

function [S_map, psi_map, q1, q2] = nematic_order_map_m(phi_rad, weight, mask, sigma_px)
% Local 2D nematic order parameter via Q-tensor averaging.
%
% Following Mirza et al. (eLife 2024, arXiv:2306.15352):
%   q1_pixel = retardance * cos(2*phi) / 2
%   q2_pixel = retardance * sin(2*phi) / 2
%   q1 = <q1_pixel>_G / <mask>_G    (area-normalized Gaussian average)
%   S = 2 * sqrt(q1^2 + q2^2)

    C  = cos(2 * phi_rad);
    Sm = sin(2 * phi_rad);

    q1_pixel = weight .* C / 2;
    q2_pixel = weight .* Sm / 2;
    q1_pixel(~mask) = 0;
    q2_pixel(~mask) = 0;

    mask_float = double(mask);
    mask_avg = imgaussfilt(mask_float, sigma_px);

    q1 = imgaussfilt(q1_pixel, sigma_px) ./ max(mask_avg, 1e-10);
    q2 = imgaussfilt(q2_pixel, sigma_px) ./ max(mask_avg, 1e-10);

    S_map = 2 * sqrt(q1.^2 + q2.^2);
    psi_map = 0.5 * atan2(q2, q1);

    S_map(~mask) = NaN;
    psi_map(~mask) = NaN;
    q1(~mask) = NaN;
    q2(~mask) = NaN;
end


function [S, psi, q1, q2] = nematic_order_global_m(phi_rad, weight, mask)
% Whole-mask scalar order parameter via Q-tensor averaging.
%
% q1 = mean(ret * cos(2*phi)) / 2
% q2 = mean(ret * sin(2*phi)) / 2
% S  = 2 * sqrt(q1^2 + q2^2)

    N = sum(mask(:));
    if N == 0
        S = NaN; psi = NaN; q1 = NaN; q2 = NaN;
        return
    end

    w  = weight(mask);
    C  = cos(2 * phi_rad(mask));
    Sm = sin(2 * phi_rad(mask));

    q1 = mean(w .* C) / 2;
    q2 = mean(w .* Sm) / 2;

    S   = 2 * sqrt(q1^2 + q2^2);
    psi = 0.5 * atan2(q2, q1);
end


function alpha = tangential_radial_alignment_m(phi_rad, mask, xc, yc)
% Angle between director and local boundary normal.
%
% Uses gradient of Euclidean distance transform for local normals.
% alpha = 0    -> radial
% alpha = pi/2 -> tangential

    [H, W] = size(mask);

    perim = mask & ~imerode(mask, strel('disk', 1));
    D = bwdist(perim);

    [Gx, Gy] = gradient(D);

    normal_angle = atan2(Gy, Gx);

    [xx, yy] = meshgrid(1:W, 1:H);
    to_center_x = xc - xx;
    to_center_y = yc - yy;
    dot_prod = Gx .* to_center_x + Gy .* to_center_y;
    flip = dot_prod < 0;
    normal_angle(flip) = normal_angle(flip) + pi;

    diff_angle = phi_rad - normal_angle;
    diff_angle = mod(diff_angle, pi);
    alpha = min(diff_angle, pi - diff_angle);

    alpha(~mask) = NaN;
end


function [S_theta, psi_theta, q1_theta, q2_theta] = ...
    cortex_nematic_profile_m(phi_rad, weight, xb, yb, xc, yc, n_theta_bins, inset_px)
% Theta-binned nematic order S(theta) around the cortex.

    [H, W] = size(phi_rad);

    % Shrink boundary inward
    dx = xb - xc;
    dy = yb - yc;
    dist = sqrt(dx.^2 + dy.^2);
    shrink = max(dist - inset_px, 1) ./ max(dist, 1e-10);
    xb_s = xc + dx .* shrink;
    yb_s = yc + dy .* shrink;

    % Angular position of each boundary point
    th = atan2(yb_s - yc, xb_s - xc);
    th(th < 0) = th(th < 0) + 2*pi;

    theta_edges = linspace(0, 2*pi, n_theta_bins + 1);

    % Interpolate fields at boundary points
    C  = cos(2 * phi_rad);
    Sm = sin(2 * phi_rad);

    [cols_grid, rows_grid] = meshgrid(1:W, 1:H);
    cvals = interp2(cols_grid, rows_grid, C,  xb_s, yb_s, 'linear', 0);
    svals = interp2(cols_grid, rows_grid, Sm, xb_s, yb_s, 'linear', 0);
    wvals = interp2(cols_grid, rows_grid, weight, xb_s, yb_s, 'linear', 0);

    bins = discretize(th, theta_edges);
    bins(isnan(bins)) = 1;

    S_theta   = nan(n_theta_bins, 1);
    psi_theta = nan(n_theta_bins, 1);
    q1_theta  = nan(n_theta_bins, 1);
    q2_theta  = nan(n_theta_bins, 1);

    for b = 1:n_theta_bins
        m = (bins == b);
        n_pts = sum(m);
        if n_pts == 0, continue; end

        q1_theta(b) = mean(wvals(m) .* cvals(m)) / 2;
        q2_theta(b) = mean(wvals(m) .* svals(m)) / 2;
        S_theta(b)  = 2 * sqrt(q1_theta(b)^2 + q2_theta(b)^2);
        psi_theta(b) = 0.5 * atan2(q2_theta(b), q1_theta(b));
    end
end


function mask = segment_from_states_m(base_dir, state_patterns, sigmaBlur, closeRadius, minArea)
% Segment oocyte from sum of 4 state images (inverted Otsu).

    Isum = 0;
    for s = 1:4
        files = dir(fullfile(base_dir, state_patterns{s}));
        files = files(~[files.isdir]);
        img = double(imread(fullfile(files(1).folder, files(1).name)));
        if size(img,3) > 1, img = mean(img,3); end
        Isum = Isum + img;
    end

    mask = segment_core_m(Isum, sigmaBlur, closeRadius, minArea, true);
end


function mask = segment_from_image_m(img, sigmaBlur, closeRadius, minArea)
% Segment oocyte from a single image. Tries both polarities.

    BW_bright = segment_core_m(img, sigmaBlur, closeRadius, minArea, false);
    BW_dark   = segment_core_m(img, sigmaBlur, closeRadius, minArea, true);

    if blob_score_m(BW_dark, size(img)) > blob_score_m(BW_bright, size(img))
        mask = BW_dark;
    else
        mask = BW_bright;
    end
end


function BW = segment_core_m(img, sigmaBlur, closeRadius, minArea, invert)
% Core segmentation: blur -> Otsu -> cleanup -> largest component.

    I_blur = imgaussfilt(double(img), sigmaBlur);
    I_norm = I_blur / max(I_blur(:) + 1e-10);
    thresh = graythresh(I_norm);

    if invert
        BW = I_norm < thresh;
    else
        BW = I_norm > thresh;
    end

    se = strel('disk', closeRadius);
    BW = imclose(BW, se);
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, minArea);

    CC = bwconncomp(BW);
    if CC.NumObjects >= 1
        numPixels = cellfun(@numel, CC.PixelIdxList);
        [~, biggest] = max(numPixels);
        BW = false(size(BW));
        BW(CC.PixelIdxList{biggest}) = true;
    end
end


function sc = blob_score_m(BW, imsize)
% Score a mask: prefer compact, centred blobs.

    if ~any(BW(:))
        sc = -1;
        return
    end
    H = imsize(1); W = imsize(2);
    [ys, xs] = find(BW);
    cx = mean(xs); cy = mean(ys);
    dist_to_center = sqrt((cx - W/2)^2 + (cy - H/2)^2) / max(H, W);
    area_frac = sum(BW(:)) / numel(BW);
    sc = area_frac - dist_to_center;
end


function [xb, yb, xc, yc, R_fit] = find_boundary_and_center_m(BW)
% Extract boundary, fit circle, return boundary coords and center.

    B = bwboundaries(BW, 'noholes');
    if isempty(B)
        xb = []; yb = []; xc = 0; yc = 0; R_fit = 0;
        return
    end

    % Find longest boundary
    lens = cellfun(@(b) size(b,1), B);
    [~, idx] = max(lens);
    bnd = B{idx};

    yb = bnd(:,1);
    xb = bnd(:,2);

    % Circle fit
    [xc, yc, R_fit] = circfit(xb, yb);
end


function save_frame_plots(phi_deg, phi_rad, weight, mask, raw_img, ...
                          xb, yb, xc, yc, R_fit, ...
                          S_map, psi_map, S_g, psi_g, ...
                          alpha_map, S_theta, psi_theta, theta_centers_deg, ...
                          q1_map, q2_map, q1_g, q2_g, q1_theta, q2_theta, ...
                          frame_dir, px_per_um, sigma_um, director_spacing)
% Generate and save analysis plots for a single frame.

    [H, W] = size(phi_deg);
    raw_img = double(raw_img);

    % 1. S heatmap
    figure('Visible','off','Position',[100 100 800 600]);
    S_disp = S_map;
    S_disp(~mask) = 0;
    imagesc(S_disp); axis image; axis off;
    colormap hot; colorbar;
    title(sprintf('S map  (\\sigma=%.1f \\mum,  S_{global}=%.4f)', sigma_um, S_g));
    saveas(gcf, fullfile(frame_dir, 'S_heatmap.png'));
    close;

    % 2. Director field on S map
    figure('Visible','off','Position',[100 100 800 600]);
    S_disp = S_map; S_disp(~mask) = 0;
    imagesc(S_disp); axis image; hold on; colormap hot;
    % Overlay directors
    step = director_spacing;
    [xx, yy] = meshgrid(step:step:W, step:step:H);
    xx = xx(:); yy = yy(:);
    valid = mask(sub2ind([H W], min(max(yy,1),H), min(max(xx,1),W)));
    xx = xx(valid); yy = yy(valid);
    for i = 1:numel(xx)
        p = psi_map(yy(i), xx(i));
        if isnan(p), continue; end
        len = step * 0.4;
        dx = len * cos(p);
        dy = len * sin(p);
        plot([xx(i)-dx, xx(i)+dx], [yy(i)-dy, yy(i)+dy], 'w-', 'LineWidth', 0.8);
    end
    axis off;
    title('Director field on S map');
    saveas(gcf, fullfile(frame_dir, 'director_on_S.png'));
    close;

    % 3. Tangential/radial alignment
    figure('Visible','off','Position',[100 100 800 600]);
    alpha_disp = rad2deg(alpha_map);
    alpha_disp(~mask) = NaN;
    imagesc(alpha_disp, [0 90]); axis image; axis off;
    colormap(cool); c = colorbar;
    c.Label.String = 'Alignment angle (deg)';
    title('Tangential (90) / Radial (0) alignment');
    saveas(gcf, fullfile(frame_dir, 'tangential_radial.png'));
    close;

    % 4. Cortex nematic profile
    figure('Visible','off','Position',[100 100 800 400]);
    subplot(2,1,1);
    plot(theta_centers_deg, S_theta, 'b-', 'LineWidth', 1.2);
    xlabel('Angle (deg)'); ylabel('S(\theta)');
    title('Cortex nematic order profile');
    xlim([0 360]);

    subplot(2,1,2);
    plot(theta_centers_deg, rad2deg(psi_theta), 'r-', 'LineWidth', 1.2);
    xlabel('Angle (deg)'); ylabel('\psi (deg)');
    title('Director angle along cortex');
    xlim([0 360]);
    saveas(gcf, fullfile(frame_dir, 'cortex_nematic_profile.png'));
    close;

    % Save per-frame results
    save(fullfile(frame_dir, 'nematic_results.mat'), ...
         'S_map', 'psi_map', 'q1_map', 'q2_map', ...
         'alpha_map', 'S_theta', 'psi_theta', 'q1_theta', 'q2_theta', ...
         'theta_centers_deg', 'S_g', 'psi_g', 'q1_g', 'q2_g');

    fprintf('  Saved frame plots to %s\n', frame_dir);
end

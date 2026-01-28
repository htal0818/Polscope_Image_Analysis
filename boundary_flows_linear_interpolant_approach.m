% Measure tangential surface flows from PIV data using LINEAR INTERPOLANT approach
%
% This script implements the simplified analytical tangent methodology:
% - Boundary points and tangents computed analytically from polynomial fit
% - PIV velocities sampled directly at boundary using griddedInterpolant
% - No upsampling from 500 to 5000 points (eliminates linear interp artifacts)
% - Tangent vectors computed analytically from polar curve derivatives
%
% Key differences from boundary_flows.m:
% 1. Tangents computed analytically: t = dr/dtheta (polar curve formula)
% 2. PIV sampled at boundary points via griddedInterpolant (linear interp of discrete PIV)
% 3. Optional offset inward along normal for cortical sampling
% 4. Direct 1:1 mapping: 500 boundary points -> 500 tangential velocities
%
% The only linear interpolation occurs where appropriate: sampling the discrete
% PIV grid at boundary locations. Geometry (boundary + tangents) stays analytical.

clear all; close all; clc

%% ========================== USER INPUTS ==================================
% --- PolScope image folders / patterns ---
base_dir = '/Users/hridaytalreja/Desktop/Jan_data_2026/jan_20_2026_FSW_and_eggs_50msexp_15sint_20x_50nmceiling/eggs/SMS_2026_0120_1518_1/Pos0/';

% If you have 4 states per frame (best for boundary):
s1 = fullfile(base_dir,'*State1*');
s2 = fullfile(base_dir,'*State2*');
s3 = fullfile(base_dir,'*State3*');
s4 = fullfile(base_dir,'*State4*');

% If you only have State1, set useFourStates=false and fill s1 accordingly.
useFourStates = true;

d1 = dir(s1);
if useFourStates
    d2 = dir(s2); d3 = dir(s3); d4 = dir(s4);
end

% --- PIVlab output .mat file ---
pivMatFile = '/Users/hridaytalreja/Desktop/Jan_data_2026/jan_20_2026_FSW_and_eggs_50msexp_15sint_20x_50nmceiling/eggs/SMS_2026_0120_1518_1/Pos0/Jan21_2026_PIV/PIVlab_output.mat';

% --- Optional crop (match your workflows) ---
doCrop = true;
cropRect = [500 500 2000 2000];  % [x y w h] in pixels

% --- Units / timing ---
px_per_um = 6.25 / 2;     % px/um
dt_sec    = 15;           % sec/frame

% --- PIV velocity units ---
pivVelUnit = 'px_per_frame';  % 'px_per_frame' or 'px_per_sec'

% --- LINEAR INTERPOLANT APPROACH PARAMETERS ---
% Offset from boundary (pixels) - sample slightly inside for cortical flow
% Set to 0 to sample exactly at boundary
normalOffsetPx = 3;  % sample 3 pixels inward along normal

% --- Curvature-based boundary parameters ---
sigmaBlur       = 1.0;           % pre-blur (pixels)
threshFrac      = 0.85;          % mask = I < threshFrac*mean2(I)
se              = strel('diamond',5);
polyOrder       = 50;            % polyfit order for r(theta)
nBoundary       = 500;           % samples along boundary (NO upsampling needed now)

% --- Theta binning for output compatibility ---
% Note: With analytical approach, we have 500 direct values at known theta
% Binning is optional - for compatibility with downstream analysis
nThetaBins = 101;
thetaEdges = linspace(0,2*pi,nThetaBins);
thetaCenters = thetaEdges;
thetaBinEdges = linspace(0,2*pi,nThetaBins+1);

% --- Quality control thresholds ---
minAreaFrac  = 0.05;
maxEccentric = 0.95;
minSolidity  = 0.85;

% --- Smart caching ---
useMaskCaching = true;
cacheIntensityThreshold = 0.02;
cacheForceRecalcEveryN  = 25;
cacheUseHintCenter      = true;

% --- Outputs ---
outDir = fullfile(base_dir, 'tangential_kymo_linear_interp');
if ~exist(outDir,'dir'); mkdir(outDir); end

outDir_kymographs = fullfile(outDir, 'kymographs');
outDir_quiver = fullfile(outDir, 'quiver_overlays');
outDir_snapshots = fullfile(outDir, 'snapshots');
outDir_qc = fullfile(outDir, 'qc');

if ~exist(outDir_kymographs,'dir'); mkdir(outDir_kymographs); end
if ~exist(outDir_quiver,'dir'); mkdir(outDir_quiver); end
if ~exist(outDir_snapshots,'dir'); mkdir(outDir_snapshots); end
if ~exist(outDir_qc,'dir'); mkdir(outDir_qc); end

saveEveryN_QC = 20;

% --- Visualization settings ---
makeQuiverOverlays = true;
quiverEveryNFrames = 10;
quiverSubsample = 3;
quiverScale = 1.5;

makeEnhancedKymographs = true;
makeSnapshotPlots = true;
snapshotEveryNFrames = 2;

%% ========================== LOAD PIV DATA =================================
S = load(pivMatFile);
[Xc, Yc, Uc, Vc] = pickPIVFields(S);

nFramesPIV = numel(Uc);
nFramesImg = numel(d1);
nFrames = min(nFramesPIV, nFramesImg);
fprintf('Frames: images=%d, piv=%d, using=%d\n', nFramesImg, nFramesPIV, nFrames);

%% ========================== PREALLOCATE ===================================
% Raw tangential velocity at 500 boundary points (no binning)
Vtheta_raw = nan(nFrames, nBoundary);       % tangential velocity at each boundary point
Theta_raw = nan(nFrames, nBoundary);        % theta values for each boundary point

% Binned output (for compatibility with original approach)
Vtheta_kymo = nan(nFrames, nThetaBins);
Npts_kymo   = zeros(nFrames, nThetaBins);

centroidXY  = nan(nFrames,2);
areaMask    = nan(nFrames,1);
qcFlag      = false(nFrames,1);
RADIUS_OF_CURVATURE = nan(nFrames, numel(thetaCenters)-1);

% Vorticity storage
Vorticity_kymo = nan(nFrames, nThetaBins);
Vorticity_field = cell(nFrames, 1);
MeanVorticity = nan(nFrames, 1);

% Cache variables
cache_prevBW = [];
cache_prevCenter = [];
cache_prevMeanInt = [];
cache_prevPoly = [];
cache_prevRADIUS = [];
cache_prevTangents = [];  % NEW: cache analytical tangents
cache_prevTheta = [];     % NEW: cache theta values
cacheHitCount = 0;
cacheRecalcCount = 0;

% Visualization data storage
visPolySeq = cell(nFrames,1);
visImageSeq = cell(nFrames,1);

%% ========================== MAIN LOOP =====================================
for fr = 1:nFrames

    %% ----- Load PolScope intensity image -----
    a1 = double(imread(fullfile(d1(fr).folder, d1(fr).name)));
    if useFourStates
        a2 = double(imread(fullfile(d2(fr).folder, d2(fr).name)));
        a3 = double(imread(fullfile(d3(fr).folder, d3(fr).name)));
        a4 = double(imread(fullfile(d4(fr).folder, d4(fr).name)));
        Iraw = a1 + a2 + a3 + a4;
    else
        Iraw = a1;
    end

    if doCrop
        Iraw = imcrop(Iraw, cropRect);
    end

    I = mat2gray(Iraw);

    %% =========================================================================
    %% SMART CACHING BLOCK
    %% =========================================================================
    needsRecalc = true;
    cacheReasonStr = 'first frame';

    if useMaskCaching && ~isempty(cache_prevBW)
        if mod(fr, cacheForceRecalcEveryN) == 1
            needsRecalc = true;
            cacheReasonStr = sprintf('forced recalc (every %d frames)', cacheForceRecalcEveryN);
        else
            meanIntCurrent = mean(I(cache_prevBW), 'omitnan');
            intensityChange = abs(meanIntCurrent - cache_prevMeanInt) / (cache_prevMeanInt + eps);

            if intensityChange < cacheIntensityThreshold
                needsRecalc = false;
                cacheReasonStr = sprintf('cache hit (delta=%.3f%%)', intensityChange*100);
            else
                needsRecalc = true;
                cacheReasonStr = sprintf('intensity changed %.3f%%', intensityChange*100);
            end
        end
    end

    if needsRecalc || ~useMaskCaching
        if ~useMaskCaching
            cacheReasonStr = 'caching disabled';
        end

        centerHint = [];
        if useMaskCaching && cacheUseHintCenter && ~isempty(cache_prevCenter)
            centerHint = cache_prevCenter;
        end

        % NEW: Get boundary with analytical tangents
        [BW, stats, poly, RADIUS, xc, yc, tx, ty, theta_boundary] = ...
            make_oocyte_mask_with_tangents(I, sigmaBlur, threshFrac, se, ...
            polyOrder, nBoundary, thetaCenters, minAreaFrac, maxEccentric, ...
            minSolidity, centerHint);

        if isempty(stats)
            fprintf('Frame %d: mask failed QC (%s).\n', fr, cacheReasonStr);
            continue;
        end

        % Update cache
        if useMaskCaching
            cache_prevBW = BW;
            cache_prevCenter = [xc, yc];
            cache_prevMeanInt = mean(I(BW), 'omitnan');
            cache_prevPoly = poly;
            cache_prevRADIUS = RADIUS;
            cache_prevTangents = [tx(:), ty(:)];  % NEW
            cache_prevTheta = theta_boundary;      % NEW
        end

        cacheRecalcCount = cacheRecalcCount + 1;

    else
        % CACHE REUSE
        BW = cache_prevBW;
        poly = cache_prevPoly;
        RADIUS = cache_prevRADIUS;
        xc = cache_prevCenter(1);
        yc = cache_prevCenter(2);
        tx = cache_prevTangents(:,1);
        ty = cache_prevTangents(:,2);
        theta_boundary = cache_prevTheta;

        stats = regionprops(BW, 'Area', 'Eccentricity', 'Centroid', 'Solidity');
        if isempty(stats)
            fprintf('Frame %d: cached mask invalid, forcing recalc.\n', fr);
            cache_prevBW = [];
            continue;
        end
        [~, ii] = max([stats.Area]);
        stats = stats(ii);

        cacheHitCount = cacheHitCount + 1;
    end

    %% =========================================================================
    %% END SMART CACHING BLOCK
    %% =========================================================================

    if isempty(stats)
        fprintf('Frame %d: mask failed QC.\n', fr);
        continue;
    end

    qcFlag(fr)    = true;
    centroidXY(fr,:) = [xc, yc];
    areaMask(fr)  = stats.Area;
    RADIUS_OF_CURVATURE(fr,:) = RADIUS;

    %% ----- Get PIV frame -----
    X = Xc{fr}; Y = Yc{fr}; U = Uc{fr}; V = Vc{fr};
    if doCrop
        X = X - cropRect(1);
        Y = Y - cropRect(2);
    end

    %% ----- Compute tangential velocity using LINEAR INTERPOLANT approach -----
    [vT_raw, theta_raw, vBins, nBinsCount, dbgFlow] = tangential_linear_interpolant( ...
        X, Y, U, V, poly, tx, ty, theta_boundary, [xc, yc], ...
        normalOffsetPx, nThetaBins, thetaBinEdges, ...
        px_per_um, dt_sec, pivVelUnit);

    % Store raw (500 points) and binned outputs
    Vtheta_raw(fr,:) = vT_raw;
    Theta_raw(fr,:) = theta_raw;
    Vtheta_kymo(fr,:) = vBins;
    Npts_kymo(fr,:)   = nBinsCount;

    %% ----- Compute vorticity -----
    [vortBins, vortField, meanVort, ~] = compute_vorticity_from_PIV( ...
        X, Y, U, V, BW, [xc, yc], ...
        1, 20, nThetaBins, thetaBinEdges, ...  % band parameters for vorticity
        px_per_um, dt_sec, pivVelUnit);

    Vorticity_kymo(fr, :) = vortBins;
    Vorticity_field{fr} = vortField;
    MeanVorticity(fr) = meanVort;

    %% ----- Store visualization data -----
    visPolySeq{fr} = poly;
    visImageSeq{fr} = I;

    %% ----- QC overlays -----
    if mod(fr, saveEveryN_QC) == 1
        qcFig = figure('Visible','off'); imshow(I,[]); hold on;
        plot(poly(:,1), poly(:,2), 'c-', 'LineWidth', 2);
        scatter(dbgFlow.sampleX, dbgFlow.sampleY, 8, 'r', 'filled');
        % Show tangent vectors at a few points
        subsamp = 1:20:nBoundary;
        quiver(poly(subsamp,1), poly(subsamp,2), tx(subsamp)*20, ty(subsamp)*20, 0, 'g', 'LineWidth', 1);
        title(sprintf('Frame %d: Boundary + Sample Points + Tangents', fr));
        legend({'Boundary', 'Sample pts', 'Tangents'}, 'Location', 'best');
        exportgraphics(qcFig, fullfile(outDir_qc, sprintf('QC_linear_interp_fr%04d.png', fr)), 'Resolution', 250);
        close(qcFig);
    end

    if mod(fr,50)==0
        fprintf('Processed frame %d/%d (%s)\n', fr, nFrames, cacheReasonStr);
    end
end

%% ========================== CACHE DIAGNOSTICS ================================
if useMaskCaching
    totalProcessed = cacheHitCount + cacheRecalcCount;
    cacheHitRate = 100 * cacheHitCount / max(1, totalProcessed);
    fprintf('\n=== SMART CACHING PERFORMANCE ===\n');
    fprintf('Cache hits:        %d/%d (%.1f%%)\n', cacheHitCount, totalProcessed, cacheHitRate);
    fprintf('Recalculations:    %d/%d (%.1f%%)\n', cacheRecalcCount, totalProcessed, 100-cacheHitRate);
    fprintf('=================================\n\n');
else
    fprintf('\nSmart caching was DISABLED\n\n');
end

%% ========================== VISUALIZATIONS ===================================
time_min = (0:nFrames-1) * (dt_sec/60);

fprintf('Generating visualizations...\n');

% --- Convert to nm/s for display (1 um/s = 1000 nm/s) ---
Vtheta_raw_nm = Vtheta_raw * 1000;
Vtheta_kymo_nm = Vtheta_kymo * 1000;

% --- 1) RAW TANGENTIAL VELOCITY (500 points, no binning) ---
fig0 = figure('Visible', 'off', 'Position', [100 100 1200 600]);
subplot(1,2,1);
imagesc(1:nBoundary, time_min, Vtheta_raw_nm);
axis tight;
xlabel('Boundary Point Index', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Raw Tangential Velocity (500 pts, analytical tangents)', 'FontSize', 12);
colormap(gca, redblue(256));
vmax = max(abs(Vtheta_raw_nm(:)), [], 'omitnan');
if ~isnan(vmax) && vmax > 0; clim([-vmax, vmax]); end
cb = colorbar; ylabel(cb, 'v_\theta (nm/s)', 'FontSize', 11);

subplot(1,2,2);
imagesc(1:nThetaBins, time_min, Vtheta_kymo_nm);
axis tight;
xlabel('Angular Bin', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Binned Tangential Velocity (for comparison)', 'FontSize', 12);
colormap(gca, redblue(256));
if ~isnan(vmax) && vmax > 0; clim([-vmax, vmax]); end
cb = colorbar; ylabel(cb, 'v_\theta (nm/s)', 'FontSize', 11);

sgtitle('LINEAR INTERPOLANT APPROACH: Raw vs Binned', 'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(fig0, fullfile(outDir_kymographs,'kymograph_raw_vs_binned.png'), 'Resolution', 300);
close(fig0);

% --- 2) BASIC SIGNED KYMOGRAPH ---
fig1 = figure('Visible', 'off', 'Position', [100 100 1000 600]);
imagesc(1:nThetaBins, time_min, Vtheta_kymo_nm);
axis tight;
xlabel('Angular Bin Number', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Tangential Cortical Flow v_\theta(\theta,t) - Linear Interpolant Approach', 'FontSize', 14);
colormap(gca, 'parula');
cb = colorbar;
ylabel(cb, 'v_\theta (nm/s)', 'FontSize', 11);
set(gca, 'FontSize', 11);
exportgraphics(fig1, fullfile(outDir_kymographs,'kymograph_vtheta_signed.png'), 'Resolution', 300);
close(fig1);

if makeEnhancedKymographs
    % --- 3) MAGNITUDE KYMOGRAPH ---
    fig2 = figure('Visible', 'off', 'Position', [100 100 1000 600]);
    imagesc(1:nThetaBins, time_min, abs(Vtheta_kymo_nm));
    axis tight;
    xlabel('Angular Bin Number', 'FontSize', 12);
    ylabel('Time (min)', 'FontSize', 12);
    title('Tangential Flow Magnitude |v_\theta| - Linear Interpolant', 'FontSize', 14);
    colormap(gca, 'hot');
    cb = colorbar;
    ylabel(cb, '|v_\theta| (nm/s)', 'FontSize', 11);
    exportgraphics(fig2, fullfile(outDir_kymographs,'kymograph_vtheta_magnitude.png'), 'Resolution', 300);
    close(fig2);

    % --- 4) DIRECTIONAL KYMOGRAPH ---
    fig3 = figure('Visible', 'off', 'Position', [100 100 1000 600]);
    imagesc(1:nThetaBins, time_min, Vtheta_kymo_nm);
    axis tight;
    xlabel('Angular Bin Number', 'FontSize', 12);
    ylabel('Time (min)', 'FontSize', 12);
    title('Tangential Flow Directionality - Linear Interpolant', 'FontSize', 14);
    colormap(gca, redblue(256));
    vmax_dir = max(abs(Vtheta_kymo_nm(:)), [], 'omitnan');
    if ~isnan(vmax_dir) && vmax_dir > 0; clim([-vmax_dir, vmax_dir]); end
    cb = colorbar;
    ylabel(cb, 'v_\theta (nm/s): red=CCW, blue=CW', 'FontSize', 10);
    exportgraphics(fig3, fullfile(outDir_kymographs,'kymograph_vtheta_directional.png'), 'Resolution', 300);
    close(fig3);

    fprintf('  - Saved enhanced kymographs\n');
end

% --- VORTICITY KYMOGRAPHS ---
fprintf('  - Generating vorticity kymographs...\n');

fig4a = figure('Visible', 'off', 'Position', [100 100 1000 600]);
imagesc(1:nThetaBins, time_min, Vorticity_kymo);
axis tight;
xlabel('Angular Bin Number', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Vorticity \omega(\theta,t) - Fluid Rotation (1/s)', 'FontSize', 14);
colormap(gca, redblue(256));
vort_max = max(abs(Vorticity_kymo(:)), [], 'omitnan');
if ~isnan(vort_max) && vort_max > 0; clim([-vort_max, vort_max]); end
cb = colorbar;
ylabel(cb, '\omega (1/s): red=CCW, blue=CW', 'FontSize', 10);
exportgraphics(fig4a, fullfile(outDir_kymographs,'kymograph_vorticity_signed.png'), 'Resolution', 300);
close(fig4a);

fig4b = figure('Visible', 'off', 'Position', [100 100 1000 600]);
imagesc(1:nThetaBins, time_min, abs(Vorticity_kymo));
axis tight;
xlabel('Angular Bin Number', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Vorticity Magnitude |\omega| (1/s)', 'FontSize', 14);
colormap(gca, 'hot');
cb = colorbar;
ylabel(cb, '|\omega| (1/s)', 'FontSize', 11);
exportgraphics(fig4b, fullfile(outDir_kymographs,'kymograph_vorticity_magnitude.png'), 'Resolution', 300);
close(fig4b);

fig4c = figure('Visible', 'off', 'Position', [100 100 800 400]);
plot(time_min, MeanVorticity, 'b-', 'LineWidth', 1.5);
hold on; yline(0, 'k--', 'LineWidth', 1);
xlabel('Time (min)', 'FontSize', 12);
ylabel('Mean Vorticity \omega (1/s)', 'FontSize', 12);
title('Global Cortical Vorticity vs Time', 'FontSize', 14);
grid on;
exportgraphics(fig4c, fullfile(outDir_kymographs,'vorticity_timeseries.png'), 'Resolution', 300);
close(fig4c);

fprintf('  - Saved vorticity kymographs\n');

% --- QUIVER OVERLAYS ---
if makeQuiverOverlays
    fprintf('  - Generating quiver overlays...\n');

    allVtheta_nm = Vtheta_kymo_nm(qcFlag, :);
    vabs_max_nm = max(abs(allVtheta_nm(:)), [], 'omitnan');
    clim_quiver = [-vabs_max_nm, vabs_max_nm];
    cmap_diverge = redblue(256);

    for fr = 1:quiverEveryNFrames:nFrames
        if ~qcFlag(fr), continue; end

        I = visImageSeq{fr};
        poly = visPolySeq{fr};
        if isempty(I) || isempty(poly), continue; end

        xc = centroidXY(fr, 1);
        yc = centroidXY(fr, 2);
        vtheta_nm = Vtheta_kymo_nm(fr, :);

        figQ = figure('Visible', 'off', 'Position', [100 100 900 800]);
        imshow(I, []); hold on; axis on;
        plot(poly(:,1), poly(:,2), 'c-', 'LineWidth', 1);

        thetaSubsample = 1:quiverSubsample:nThetaBins;
        nQuiver = numel(thetaSubsample);

        xq = zeros(nQuiver, 1);
        yq = zeros(nQuiver, 1);
        uq = zeros(nQuiver, 1);
        vq = zeros(nQuiver, 1);
        vq_color = zeros(nQuiver, 1);

        for k = 1:nQuiver
            idx = thetaSubsample(k);
            theta_k = thetaCenters(idx);
            vtheta_k_nm = vtheta_nm(idx);

            if isnan(vtheta_k_nm)
                vq_color(k) = NaN;
                continue;
            end

            R_mean = mean(RADIUS_OF_CURVATURE(fr, :), 'omitnan');
            if isnan(R_mean), R_mean = 100; end

            xq(k) = xc + R_mean * cos(theta_k);
            yq(k) = yc + R_mean * sin(theta_k);

            % Tangent vector (perpendicular to radial)
            tx_k = -sin(theta_k);
            ty_k = cos(theta_k);

            arrow_scale = quiverScale * abs(vtheta_k_nm) / (px_per_um * 1000);  % nm/s to display
            if vtheta_k_nm < 0
                tx_k = -tx_k;
                ty_k = -ty_k;
            end

            uq(k) = arrow_scale * tx_k;
            vq(k) = arrow_scale * ty_k;
            vq_color(k) = vtheta_k_nm;
        end

        valid = ~isnan(vq_color) & (xq ~= 0 | yq ~= 0);
        xq = xq(valid); yq = yq(valid);
        uq = uq(valid); vq = vq(valid);
        vq_color = vq_color(valid);

        nColors = size(cmap_diverge, 1);
        colorIdx = round((vq_color - clim_quiver(1)) / (clim_quiver(2) - clim_quiver(1)) * (nColors - 1)) + 1;
        colorIdx = max(1, min(nColors, colorIdx));

        for k = 1:numel(xq)
            arrowColor = cmap_diverge(colorIdx(k), :);
            quiver(xq(k), yq(k), uq(k), vq(k), 0, 'Color', arrowColor, 'LineWidth', 1.5, 'MaxHeadSize', 2);
        end

        hScatter = scatter(xq, yq, 1, vq_color, 'filled');
        set(hScatter, 'MarkerFaceAlpha', 0);
        colormap(gca, cmap_diverge);
        caxis(clim_quiver);

        cb = colorbar;
        ylabel(cb, 'v_\theta (nm/s)', 'FontSize', 10, 'Color', 'w');
        set(cb, 'Color', 'w');

        title(sprintf('Frame %d: Tangential Flow (Linear Interp)', fr), 'FontSize', 12, 'Color', 'w');

        exportgraphics(figQ, fullfile(outDir_quiver, sprintf('quiver_fr%04d.png', fr)), 'Resolution', 200);
        close(figQ);
    end

    fprintf('  - Saved quiver overlays\n');
end

% --- SNAPSHOT PLOTS ---
if makeSnapshotPlots
    fprintf('  - Generating snapshot plots...\n');

    % Compute global axis limits (in nm/s)
    vtheta_all_nm = Vtheta_kymo_nm(qcFlag, :);
    vtheta_absmax_nm = max(abs(vtheta_all_nm(:)), [], 'omitnan');
    if isnan(vtheta_absmax_nm) || vtheta_absmax_nm == 0, vtheta_absmax_nm = 1; end
    ylim_vtheta = [-vtheta_absmax_nm, vtheta_absmax_nm] * 1.1;

    vort_all = Vorticity_kymo(qcFlag, :);
    vort_absmax = max(abs(vort_all(:)), [], 'omitnan');
    if isnan(vort_absmax) || vort_absmax == 0, vort_absmax = 1; end
    ylim_vort = [-vort_absmax, vort_absmax] * 1.1;

    % Curvature in um (convert from px)
    curv_all_um = RADIUS_OF_CURVATURE(qcFlag, :) / px_per_um;
    curv_min_um = min(curv_all_um(:), [], 'omitnan');
    curv_max_um = max(curv_all_um(:), [], 'omitnan');
    if isnan(curv_min_um), curv_min_um = 0; end
    if isnan(curv_max_um) || curv_max_um == 0, curv_max_um = 1; end
    ylim_curv_um = [curv_min_um * 0.9, curv_max_um * 1.1];

    for fr = 1:snapshotEveryNFrames:nFrames
        if ~qcFlag(fr), continue; end

        I = visImageSeq{fr};
        poly = visPolySeq{fr};
        if isempty(I) || isempty(poly), continue; end

        xc = centroidXY(fr, 1);
        yc = centroidXY(fr, 2);
        vtheta_nm = Vtheta_kymo_nm(fr, :);
        radius_curv_um = RADIUS_OF_CURVATURE(fr, :) / px_per_um;  % Convert to um
        vort = Vorticity_kymo(fr, :);

        % Compute mask radius in um from area
        mask_radius_um = sqrt(areaMask(fr) / pi) / px_per_um;

        figSnap = figure('Visible', 'off', 'Position', [50 50 1400 900]);

        % Panel 1: Image with boundary
        subplot(2,3,1);
        imshow(I, []); hold on;
        plot(poly(:,1), poly(:,2), 'y-', 'LineWidth', 2);
        plot(xc, yc, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
        title(sprintf('Frame %d (t=%.2f min)', fr, time_min(fr)));

        % Panel 2: Binned tangential velocity (nm/s)
        subplot(2,3,2);
        plot(1:nThetaBins, vtheta_nm, 'b-', 'LineWidth', 1.5);
        hold on; yline(0, 'k--', 'LineWidth', 0.5);
        xlabel('Angular Bin'); ylabel('v_\theta (nm/s)');
        title('Tangential Velocity');
        grid on; xlim([1 nThetaBins]); ylim(ylim_vtheta);

        % Panel 3: Curvature (um)
        subplot(2,3,3);
        plot(1:nThetaBins-1, radius_curv_um, 'r-', 'LineWidth', 1.5);
        xlabel('Angular Bin'); ylabel('R_{curv} (\mum)');
        title('Boundary Curvature');
        grid on; xlim([1 nThetaBins-1]); ylim(ylim_curv_um);

        % Panel 4: Vorticity (1/s)
        subplot(2,3,4);
        plot(1:nThetaBins, vort, 'm-', 'LineWidth', 1.5);
        hold on; yline(0, 'k--', 'LineWidth', 0.5);
        xlabel('Angular Bin'); ylabel('\omega (1/s)');
        title('Vorticity Profile');
        grid on; xlim([1 nThetaBins]); ylim(ylim_vort);

        % Panel 5: Polar magnitude (nm/s)
        subplot(2,3,5);
        vtheta_abs_nm = abs(vtheta_nm);
        vtheta_abs_nm(isnan(vtheta_abs_nm)) = 0;
        polarplot(thetaCenters, vtheta_abs_nm, 'b-', 'LineWidth', 2);
        title('|v_\theta| Polar (nm/s)');
        rlim([0, vtheta_absmax_nm * 1.1]);

        % Panel 6: Statistics
        subplot(2,3,6); axis off;
        vtheta_valid_nm = vtheta_nm(~isnan(vtheta_nm));

        stats_text = {
            sprintf('Frame: %d', fr)
            sprintf('Time: %.2f min', time_min(fr))
            ''
            sprintf('v_\\theta mean: %.2f nm/s', mean(vtheta_valid_nm))
            sprintf('v_\\theta std:  %.2f nm/s', std(vtheta_valid_nm))
            ''
            sprintf('Mask radius: %.1f \\mum', mask_radius_um)
            sprintf('Normal offset: %d px', normalOffsetPx)
        };

        text(0.1, 0.9, stats_text, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'FixedWidth');

        sgtitle(sprintf('Linear Interpolant Analysis - Frame %d', fr), 'FontSize', 14, 'FontWeight', 'bold');

        exportgraphics(figSnap, fullfile(outDir_snapshots, sprintf('snapshot_fr%04d.png', fr)), 'Resolution', 200);
        close(figSnap);
    end

    fprintf('  - Saved snapshot plots\n');
end

fprintf('All visualizations complete!\n\n');

%% ========================== SAVE RESULTS ===================================
save(fullfile(outDir,'tangential_linear_interp_results.mat'), ...
     'Vtheta_raw', 'Theta_raw', ...  % NEW: raw 500-point data
     'Vtheta_kymo', 'Npts_kymo', 'thetaCenters', 'thetaBinEdges', 'time_min', ...
     'centroidXY', 'areaMask', 'qcFlag', 'RADIUS_OF_CURVATURE', ...
     'Vorticity_kymo', 'Vorticity_field', 'MeanVorticity', ...
     'px_per_um', 'dt_sec', 'normalOffsetPx', 'pivMatFile', 'base_dir');

fprintf('Saved outputs to: %s\n', outDir);
fprintf('\nLINEAR INTERPOLANT APPROACH outputs:\n');
fprintf('  - Vtheta_raw: Raw tangential velocity at 500 boundary points\n');
fprintf('  - Theta_raw: Theta values for each boundary point\n');
fprintf('  - Vtheta_kymo: Binned tangential velocity (for compatibility)\n');
fprintf('  - normalOffsetPx: %d pixels inward from boundary\n', normalOffsetPx);


%% ============================== FUNCTIONS =================================

function [BW, stats, poly, RADIUS, xc, yc, tx, ty, theta_boundary] = ...
    make_oocyte_mask_with_tangents(I, sigmaBlur, threshFrac, se, polyOrder, ...
    nBoundary, theta, minAreaFrac, maxEccentric, minSolidity, centerHint)
% Curvature-based boundary reconstruction WITH ANALYTICAL TANGENTS
%
% This function extends make_oocyte_mask_curvature to also output:
%   tx, ty - unit tangent vectors at each boundary point (analytical)
%   theta_boundary - theta values at each boundary point
%
% Tangents are computed analytically from the polar curve formula:
%   dx/dtheta = r'*cos(theta) - r*sin(theta)
%   dy/dtheta = r'*sin(theta) + r*cos(theta)

[H, W] = size(I);

if nargin < 11 || isempty(centerHint)
    centerHint = [W/2, H/2];
end

% Default outputs
BW = false(size(I));
stats = [];
poly = [];
RADIUS = nan(1, numel(theta)-1);
xc = centerHint(1); yc = centerHint(2);
tx = []; ty = []; theta_boundary = [];

% 1) Coarse threshold mask
Iblur = imgaussfilt(I, sigmaBlur);
BW0 = Iblur < (threshFrac * mean2(Iblur));
BW0 = imdilate(BW0, se);
BW0 = imfill(BW0, 'holes');
BW0 = imerode(BW0, se);

% Keep largest connected component
labels = bwlabel(BW0);
if max(labels(:)) == 0
    return;
end
Area = zeros(1, max(labels(:)));
for i = 1:max(labels(:))
    temp = regionprops(labels == i, 'Area');
    Area(i) = temp.Area;
end
mask = labels == find(Area == max(Area), 1);
BW0 = mask;

% 2) Edge pixels
edges = imgradient(BW0);
edges = edges > 0;

xx = []; yy = [];
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        if edges(i, j) == 1
            xx = [xx j];
            yy = [yy i];
        end
    end
end

if numel(xx) < 50
    warning('Too few edge pixels for circfit; returning empty mask.');
    return;
end

% 3) Circle fit
[~, xc, yc] = circfit(xx, yy);

% 4) Polar coordinates
nPts = numel(xx);
r = zeros(1, nPts);
angle = zeros(1, nPts);

for kk = 1:nPts
    r(kk) = norm([xx(kk) - xc, yy(kk) - yc]);
    if yy(kk) > yc
        angle(kk) = acos(dot(([xx(kk) - xc, yy(kk) - yc]) / norm([xx(kk) - xc, yy(kk) - yc]), [1 0]));
    else
        angle(kk) = 2*pi - acos(dot(([xx(kk) - xc, yy(kk) - yc]) / norm([xx(kk) - xc, yy(kk) - yc]), [1 0]));
    end
end

% 5) FIRST PASS: polynomial fit
RRrow = nan(1, nBoundary);

param = polyfit(angle, r, polyOrder);
xgrid = linspace(0, 2*pi, nBoundary);
y1 = polyval(param, xgrid);

r_theta_p = polyder(param);
r_2theta_p = polyder(r_theta_p);
r_theta = polyval(r_theta_p, xgrid);
r_2theta = polyval(r_2theta_p, xgrid);

% Curvature for middle half
for ii = (length(theta)-1)/2 - (length(theta)-1)/4 : (length(theta)-1)/2 + (length(theta)-1)/4
    sel = xgrid > theta(ii) & xgrid < theta(ii+1);
    if any(sel)
        num = ((y1(sel).^2 + r_theta(sel).^2).^(3/2));
        den = abs(y1(sel).^2 + 2*r_theta(sel).^2 - y1(sel).*r_2theta(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

RRrow(126:375) = y1(126:375);

% 6) SECOND PASS for wrap-around
[angle2, Iord] = sort(angle);
r2 = r(Iord);

angle2 = angle2 + pi;
angle2(angle2 > 2*pi) = angle2(angle2 > 2*pi) - 2*pi;

param2 = polyfit(angle2, r2, polyOrder);
x2 = linspace(0, 2*pi, nBoundary);
y2 = polyval(param2, x2);

r_theta_p2 = polyder(param2);
r_2theta_p2 = polyder(r_theta_p2);
r_theta2 = polyval(r_theta_p2, x2);
r_2theta2 = polyval(r_2theta_p2, x2);

x2 = x2 - pi;
x2(x2 < 0) = 2*pi + x2(x2 < 0);

y2 = circshift(y2, nBoundary/2);

RRrow(1:125) = y2(1:125);
RRrow(376:500) = y2(376:500);

% Curvature for remaining quarters
for ii = 1:(length(theta)-1)/2 - (length(theta)-1)/4
    sel = x2 > theta(ii) & x2 < theta(ii+1);
    if any(sel)
        num = ((y2(sel).^2 + r_theta2(sel).^2).^(3/2));
        den = abs(y2(sel).^2 + 2*r_theta2(sel).^2 - y2(sel).*r_2theta2(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

for ii = (length(theta)-1)/2 + (length(theta)-1)/4 : length(theta)-1
    sel = x2 > theta(ii) & x2 < theta(ii+1);
    if any(sel)
        num = ((y2(sel).^2 + r_theta2(sel).^2).^(3/2));
        den = abs(y2(sel).^2 + 2*r_theta2(sel).^2 - y2(sel).*r_2theta2(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

% 7) Reconstruct boundary
theta_boundary = linspace(0, 2*pi, nBoundary);
XX = RRrow .* cos(theta_boundary) + xc;
YY = RRrow .* sin(theta_boundary) + yc;

XX = min(max(XX, 1), W);
YY = min(max(YY, 1), H);

BW = poly2mask(XX, YY, H, W);
BW = imfill(BW, 'holes');

poly = [XX(:) YY(:)];

% =========================================================================
% NEW: COMPUTE ANALYTICAL TANGENTS
% =========================================================================
% For a polar curve r(theta), the tangent direction is:
%   dx/dtheta = r'(theta)*cos(theta) - r(theta)*sin(theta)
%   dy/dtheta = r'(theta)*sin(theta) + r(theta)*cos(theta)
%
% We need to use the appropriate polynomial for each region:
%   - Middle region (indices 126:375): use param (first pass)
%   - Edge regions (1:125, 376:500): use param2 (second pass)

tx = zeros(nBoundary, 1);
ty = zeros(nBoundary, 1);

% Middle region: use first pass polynomial (param)
middle_idx = 126:375;
theta_mid = theta_boundary(middle_idx);
r_mid = polyval(param, theta_mid);
r_prime_mid = polyval(r_theta_p, theta_mid);

dx_mid = r_prime_mid .* cos(theta_mid) - r_mid .* sin(theta_mid);
dy_mid = r_prime_mid .* sin(theta_mid) + r_mid .* cos(theta_mid);
mag_mid = sqrt(dx_mid.^2 + dy_mid.^2) + eps;
tx(middle_idx) = dx_mid ./ mag_mid;
ty(middle_idx) = dy_mid ./ mag_mid;

% Edge regions: use second pass polynomial (param2)
% Need to shift theta by pi for the second pass polynomial
edge_idx = [1:125, 376:500];
theta_edge = theta_boundary(edge_idx);

% Shift theta for param2 evaluation (param2 was fit with angles shifted by pi)
theta_edge_shifted = theta_edge + pi;
theta_edge_shifted(theta_edge_shifted > 2*pi) = theta_edge_shifted(theta_edge_shifted > 2*pi) - 2*pi;

r_edge = polyval(param2, theta_edge_shifted);
r_prime_edge = polyval(r_theta_p2, theta_edge_shifted);

% Compute tangent in the shifted coordinate system, then transform back
% The tangent direction in original theta is still computed using original theta
dx_edge = r_prime_edge .* cos(theta_edge) - r_edge .* sin(theta_edge);
dy_edge = r_prime_edge .* sin(theta_edge) + r_edge .* cos(theta_edge);
mag_edge = sqrt(dx_edge.^2 + dy_edge.^2) + eps;
tx(edge_idx) = dx_edge ./ mag_edge;
ty(edge_idx) = dy_edge ./ mag_edge;

% =========================================================================
% END ANALYTICAL TANGENTS
% =========================================================================

% QC
st = regionprops(BW, 'Area', 'Eccentricity', 'Centroid', 'Solidity');
if isempty(st)
    BW = false(size(BW));
    stats = [];
    poly = [];
    tx = []; ty = []; theta_boundary = [];
    return;
end
[~, ii] = max([st.Area]);
stats = st(ii);

imgArea = numel(I);
if stats.Area < minAreaFrac*imgArea || stats.Eccentricity > maxEccentric || stats.Solidity < minSolidity
    BW = false(size(BW));
    stats = [];
    poly = [];
    tx = []; ty = []; theta_boundary = [];
    return;
end

end


function [vT_raw, theta_raw, vBins, nCount, dbg] = tangential_linear_interpolant( ...
    X, Y, U, V, poly, tx, ty, theta_boundary, centroid, ...
    normalOffsetPx, nThetaBins, thetaBinEdges, ...
    px_per_um, dt_sec, pivVelUnit)
% TANGENTIAL_LINEAR_INTERPOLANT - Compute tangential velocity using griddedInterpolant
%
% This implements the simplified analytical approach:
%   1. Boundary points (xb, yb) are already computed from polynomial
%   2. Tangent vectors (tx, ty) are pre-computed analytically
%   3. PIV velocities are sampled at boundary points using griddedInterpolant
%   4. Tangential velocity = dot(velocity, tangent)
%
% INPUTS:
%   X, Y, U, V       - PIV grid and velocities
%   poly             - Boundary polygon [N x 2]
%   tx, ty           - Pre-computed analytical unit tangent vectors [N x 1]
%   theta_boundary   - Theta values at boundary points [1 x N]
%   centroid         - [xc, yc] center
%   normalOffsetPx   - Pixels to offset inward along normal (0 = on boundary)
%   nThetaBins       - Number of angular bins for binned output
%   thetaBinEdges    - Bin edges for angular binning
%   px_per_um, dt_sec, pivVelUnit - Unit conversion parameters
%
% OUTPUTS:
%   vT_raw           - Raw tangential velocity at each boundary point [1 x N]
%   theta_raw        - Theta at each boundary point [1 x N]
%   vBins            - Binned tangential velocity [1 x nThetaBins]
%   nCount           - Count per bin [1 x nThetaBins]
%   dbg              - Debug structure

% Initialize outputs
nBoundary = size(poly, 1);
vT_raw = nan(1, nBoundary);
theta_raw = theta_boundary(:)';
vBins = nan(1, nThetaBins);
nCount = zeros(1, nThetaBins);
dbg = struct('sampleX', [], 'sampleY', [], 'tx', tx, 'ty', ty);

if isempty(poly) || isempty(tx)
    return;
end

% Convert velocities to um/s
[U_um_s, V_um_s] = convertVelToUmPerSec(U, V, pivVelUnit, px_per_um, dt_sec);

% Get boundary points
xb = poly(:,1);
yb = poly(:,2);

% Compute sample points (optionally offset inward along normal)
if normalOffsetPx > 0
    % Normal is perpendicular to tangent (inward pointing)
    % For CCW boundary: normal = (-ty, tx) points inward
    nx = -ty(:);
    ny = tx(:);

    x_sample = xb + normalOffsetPx * nx;
    y_sample = yb + normalOffsetPx * ny;
else
    x_sample = xb;
    y_sample = yb;
end

dbg.sampleX = x_sample;
dbg.sampleY = y_sample;

% Create griddedInterpolant for PIV velocities
% Extract unique grid vectors from meshgrid matrices
x_vec = X(1, :);    % 1 x nCols
y_vec = Y(:, 1)';   % 1 x nRows

% Handle NaN values - griddedInterpolant doesn't like NaN in values
U_interp = U_um_s;
V_interp = V_um_s;
U_interp(isnan(U_interp)) = 0;  % or use inpaint_nans if available
V_interp(isnan(V_interp)) = 0;

% Create interpolants
% Note: griddedInterpolant expects {y_vec, x_vec} for row-column indexing
try
    F_u = griddedInterpolant({y_vec, x_vec}, U_interp, 'linear', 'none');
    F_v = griddedInterpolant({y_vec, x_vec}, V_interp, 'linear', 'none');
catch
    % Fallback to interp2 if griddedInterpolant fails
    warning('griddedInterpolant failed, using interp2 fallback');
    u_at_boundary = interp2(X, Y, U_um_s, x_sample, y_sample, 'linear');
    v_at_boundary = interp2(X, Y, V_um_s, x_sample, y_sample, 'linear');

    % Compute tangential velocity
    vT_raw = (u_at_boundary(:) .* tx(:) + v_at_boundary(:) .* ty(:))';

    % Bin for compatibility
    binIdx = discretize(theta_raw, thetaBinEdges);
    for k = 1:nThetaBins
        m = (binIdx == k) & ~isnan(vT_raw);
        nCount(k) = sum(m);
        if any(m), vBins(k) = mean(vT_raw(m), 'omitnan'); end
    end
    return;
end

% Sample velocities at boundary points
u_at_boundary = F_u(y_sample, x_sample);
v_at_boundary = F_v(y_sample, x_sample);

% Compute tangential velocity: vT = u*tx + v*ty
vT_raw = (u_at_boundary(:) .* tx(:) + v_at_boundary(:) .* ty(:))';

% Mark points outside PIV grid as NaN
outsideGrid = isnan(u_at_boundary) | isnan(v_at_boundary);
vT_raw(outsideGrid) = NaN;

% Bin for compatibility with original approach
binIdx = discretize(theta_raw, thetaBinEdges);
for k = 1:nThetaBins
    m = (binIdx == k) & ~isnan(vT_raw);
    nCount(k) = sum(m);
    if any(m), vBins(k) = mean(vT_raw(m), 'omitnan'); end
end

end


function [U_um_s, V_um_s] = convertVelToUmPerSec(U, V, pivVelUnit, px_per_um, dt_sec)
switch lower(strtrim(pivVelUnit))
    case 'px_per_frame'
        U_um_s = (U / px_per_um) / dt_sec;
        V_um_s = (V / px_per_um) / dt_sec;
    case 'px_per_sec'
        U_um_s = (U / px_per_um);
        V_um_s = (V / px_per_um);
    otherwise
        error('Unknown pivVelUnit: %s', pivVelUnit);
end
end


function [Xc, Yc, Uc, Vc] = pickPIVFields(S)
if isfield(S,'X') && isfield(S,'Y')
    Xc = S.X; Yc = S.Y;
elseif isfield(S,'x') && isfield(S,'y')
    Xc = S.x; Yc = S.y;
else
    error('Cannot find X/Y in PIV .mat');
end

if isfield(S,'U') && isfield(S,'V')
    Uc = S.U; Vc = S.V;
elseif isfield(S,'u') && isfield(S,'v')
    Uc = S.u; Vc = S.v;
elseif isfield(S,'u_original') && isfield(S,'v_original')
    Uc = S.u_original; Vc = S.v_original;
else
    error('Cannot find U/V in PIV .mat');
end

if ~iscell(Uc), Uc={Uc}; end
if ~iscell(Vc), Vc={Vc}; end
if ~iscell(Xc), Xc={Xc}; end
if ~iscell(Yc), Yc={Yc}; end

if numel(Xc)==1 && numel(Uc)>1, Xc = repmat(Xc, size(Uc)); end
if numel(Yc)==1 && numel(Uc)>1, Yc = repmat(Yc, size(Uc)); end
end


function [vortBins, vortField, meanVort, dbgVort] = compute_vorticity_from_PIV( ...
    X, Y, U, V, BW, centroid, bandOuterPx, bandInnerPx, nThetaBins, thetaBinEdges, ...
    px_per_um, dt_sec, pivVelUnit)
% Compute vorticity field from PIV data
% Vorticity: omega = dv/dx - du/dy [1/s]

vortBins = nan(1, nThetaBins);
vortField = nan(size(U));
meanVort = NaN;
dbgVort = struct('sampleX', [], 'sampleY', [], 'sampleVort', []);

[U_um_s, V_um_s] = convertVelToUmPerSec(U, V, pivVelUnit, px_per_um, dt_sec);

% Grid spacing
if size(X, 2) > 1
    dx_px = abs(X(1, 2) - X(1, 1));
else
    dx_px = 1;
end
if size(Y, 1) > 1
    dy_px = abs(Y(2, 1) - Y(1, 1));
else
    dy_px = 1;
end

dx_um = dx_px / px_per_um;
dy_um = dy_px / px_per_um;

[nRows, nCols] = size(U_um_s);

% dv/dx
dvdx = zeros(nRows, nCols);
for i = 1:nRows
    for j = 2:nCols-1
        dvdx(i, j) = (V_um_s(i, j+1) - V_um_s(i, j-1)) / (2 * dx_um);
    end
    if nCols >= 2
        dvdx(i, 1) = (V_um_s(i, 2) - V_um_s(i, 1)) / dx_um;
        dvdx(i, nCols) = (V_um_s(i, nCols) - V_um_s(i, nCols-1)) / dx_um;
    end
end

% du/dy
dudy = zeros(nRows, nCols);
for j = 1:nCols
    for i = 2:nRows-1
        dudy(i, j) = (U_um_s(i+1, j) - U_um_s(i-1, j)) / (2 * dy_um);
    end
    if nRows >= 2
        dudy(1, j) = (U_um_s(2, j) - U_um_s(1, j)) / dy_um;
        dudy(nRows, j) = (U_um_s(nRows, j) - U_um_s(nRows-1, j)) / dy_um;
    end
end

vortField = dvdx - dudy;
nanMask = isnan(U_um_s) | isnan(V_um_s);
vortField(nanMask) = NaN;

% Sample in cortical band
per = bwperim(BW);
D = bwdist(per);

Xi = clamp(round(X), 1, size(D, 2));
Yi = clamp(round(Y), 1, size(D, 1));
lin = sub2ind(size(D), Yi, Xi);

inside = BW(lin) & isfinite(vortField);
inBand = inside & (D(lin) >= bandOuterPx) & (D(lin) <= bandInnerPx);

if ~any(inBand(:))
    return;
end

xq = X(inBand);
yq = Y(inBand);
vortq = vortField(inBand);

dbgVort.sampleX = xq;
dbgVort.sampleY = yq;
dbgVort.sampleVort = vortq;

cx = centroid(1);
cy = centroid(2);
th = atan2(yq - cy, xq - cx);
th = wrapTo2Pi(th);

binIdx = discretize(th, thetaBinEdges);
for k = 1:nThetaBins
    m = (binIdx == k);
    if any(m)
        vortBins(k) = mean(vortq(m), 'omitnan');
    end
end

meanVort = mean(vortq, 'omitnan');

end


function v = clamp(v, lo, hi)
v = max(lo, min(hi, v));
end


function cmap = redblue(n)
% Diverging red-white-blue colormap
if nargin < 1, n = 256; end
r = linspace(0, 1, n)';
cmap = [r, 1-abs(2*r-1), 1-r];
end

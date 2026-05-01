% wavefront_kymograph.m
% Per-cortex-sector (depth x time) kymographs — a tiled grid of inward
% retardance profiles, one panel per angular sector around the cortex.
%
% Why this is a separate script from contour_retardance.m:
%   contour_retardance.m's distKymo angle-averages every interior pixel
%   into a single (depth x time) plot, which dilutes a wave living in one
%   cortex sector with all the quiescent sectors at the same depth — the
%   localized inward wavefront washes out. This script keeps the angle
%   resolution and presents the result as small multiples: one
%   (depth x time) kymograph per cortex sector. Wave-carrying sectors
%   light up; quiescent sectors show only the cortex band. Reading across
%   panels reveals where the wave is and how its arrival shifts around
%   the cortex over time, with no angular averaging anywhere.
%
% The map per frame is computed in one vectorized accumarray pass over
% the interior pixels using the distance transform — no per-angle ray
% loop, no per-boundary-point loop. The full (angle x depth x frame)
% cube is stored so the user can re-render with a different choice of
% panels offline without re-running segmentation.
%
% Inputs are the same image stacks consumed by contour_retardance.m
% (4-state Polscope is the default; retardance and multipage are also
% supported).
%
% Outputs (in outDir):
%   radial_kymos_grid.png/.fig             tiled grid: nPanels (depth x time) kymographs at evenly spaced cortex sectors
%   wavefront_kymograph_results.mat        full sliceCube + axes + diagnostics
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
maxDepth_um  = 40;      % how far inward from cortex to sample (microns)

% --- Tiled-grid plot ---
% Number of cortex sectors to display as separate (depth x time) panels.
% 9 -> 3x3 grid; small enough to read each panel, dense enough to see how
% wave arrival shifts panel-to-panel around the cortex.
nPanels = 9;

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

% Per-frame outputs: full (angle x depth x frame) cube + diagnostics.
sliceCube           = nan(nAngleBins, nDepth, nFrames);  % all-angles, all-depths, all-frames retardance map
avgInward_nm        = nan(nFrames, nDepth);              % per-frame angular mean (cross-check vs distKymo)
centroidXY          = nan(nFrames, 2);
fitRadius_um        = nan(nFrames, 1);

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

    % Store the full angle-resolved map; no reduction across angles. The
    % grid plot below picks K rows from this cube to render as separate
    % (depth x time) kymographs.
    sliceCube(:, :, fr) = slice;

    % Per-frame angular mean — the rotationally symmetric part the cortex
    % band lives in. Saved for cross-checking against contour_retardance's
    % distKymo, which averages over all angles and so wash-out the wave.
    avgInward_nm(fr, :) = mean(slice, 1, 'omitnan');

    if mod(fr, max(1, round(nFrames/10))) == 0
        fprintf('  Processed %d / %d frames (%.0f%%)\n', fr, nFrames, 100*fr/nFrames);
    end
end

elapsed = toc;
fprintf('Done! %.1f sec total (%.2f sec/frame)\n\n', elapsed, elapsed/nFrames);

%% ========================== PLOTS =========================================
fprintf('Generating plots...\n');

% Pick nPanels evenly spaced cortex sectors. The +1 / drop-end pattern
% avoids putting the same panel at both 0 and 2*pi.
panelIdx        = round(linspace(1, nAngleBins+1, nPanels+1));
panelIdx(end)   = [];
panelIdx        = unique(min(max(panelIdx, 1), nAngleBins));
nPanelsActual   = numel(panelIdx);
panelAngles_deg = rad2deg(angleCenters(panelIdx));

% Pull the K (depth x time) panels out of the cube.
panelData = squeeze(sliceCube(panelIdx, :, :));   % nPanels x nDepth x nFrames
panelData = permute(panelData, [3 2 1]);          % nFrames x nDepth x nPanels

% Shared color limits so panel amplitudes are comparable. Robust quantiles
% avoid having one outlier panel saturate the colormap of the others.
flat = panelData(~isnan(panelData));
if isempty(flat)
    clim_lo = 0; clim_hi = 1;
else
    clim_lo = quantile(flat, 0.02);
    clim_hi = quantile(flat, 0.98);
    if clim_hi <= clim_lo; clim_hi = clim_lo + eps; end
end

nRows = ceil(sqrt(nPanelsActual));
nCols = ceil(nPanelsActual / nRows);

fig = figure('Position', [100 100 1400 1100]);
t = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nPanelsActual
    nexttile;
    imagesc(depthAxis_um, time_min, panelData(:, :, k), [clim_lo clim_hi]);
    set(gca, 'YDir', 'normal');
    title(sprintf('\\theta = %.0f^{\\circ}', panelAngles_deg(k)));
end
xlabel(t, 'Depth from cortex (\mum)', 'FontSize', 12);
ylabel(t, 'Time (min)', 'FontSize', 12);
title(t, 'Inward retardance kymograph at each cortex sector', 'FontSize', 14);
colormap(parula);
cb = colorbar;
cb.Layout.Tile  = 'east';
cb.Label.String = 'Retardance (nm)';
cb.Label.FontSize = 11;
exportgraphics(fig, fullfile(outDir, 'radial_kymos_grid.png'), 'Resolution', 200);
savefig(fig, fullfile(outDir, 'radial_kymos_grid.fig'));
close(fig);

%% ========================== SAVE DATA =====================================
results = struct();
results.sliceCube           = sliceCube;             % full angle x depth x frame retardance map (nm)
results.avgInward_nm        = avgInward_nm;          % per-frame angular mean (cross-check vs distKymo)
results.angleCenters_rad    = angleCenters;
results.angleEdges_rad      = angleEdges;
results.depthAxis_um        = depthAxis_um;
results.depthEdges_um       = depthEdges;
results.time_sec            = time_sec;
results.time_min            = time_min;
results.centroidXY          = centroidXY;
results.fitRadius_um        = fitRadius_um;
results.nPanels             = nPanelsActual;
results.panelIdx            = panelIdx;
results.panelAngles_deg     = panelAngles_deg;
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

save(fullfile(outDir, 'wavefront_kymograph_results.mat'), '-struct', 'results');
fprintf('Saved results to %s\n', fullfile(outDir, 'wavefront_kymograph_results.mat'));

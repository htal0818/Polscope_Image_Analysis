%% av_axis_analysis.m
% Sample retardance and PIV flow along a user-defined AV (animal-vegetal)
% axis of an oocyte. Single-file, no helpers.
%
% WORKFLOW
%   1. Make sure PIV data is in the workspace (u_original, v_original, x, y),
%      same as for boundary_flows_linear_interpolant_approach.m.
%   2. Set the variables under CONFIG below (or override in the workspace).
%   3. Run the script. A figure pops up with the first retardance frame;
%      click *two* points:
%        P1 = animal pole  (boundary point near the nucleus / GV)
%        P2 = vegetal pole (opposite boundary point, or wherever you want
%                            the axis to end)
%   4. Three kymographs (position-along-AV x time) are written to outDir:
%        - retardance_kymo   (nm)
%        - flow_along_kymo   (um/s, + = toward P2)
%        - flow_perp_kymo    (um/s, + = 90 deg CCW from AV)
%      plus av_axis_results.mat and av_axis_summary.png.

%% ========================== CONFIG ==========================
if ~exist('base_dir','var') || isempty(base_dir)
    base_dir = pwd;
end
if ~exist('retardanceGlob','var') || isempty(retardanceGlob)
    retardanceGlob = '*Ret*.tif*';
end
if ~exist('px_per_um','var'); px_per_um  = 6.25 / 2; end
if ~exist('dt_sec','var');    dt_sec     = 15;       end
if ~exist('pivVelUnit','var');pivVelUnit = 'px_per_frame'; end
if ~exist('outDir','var') || isempty(outDir)
    outDir = fullfile(base_dir, 'av_axis_out');
end
if ~exist(outDir, 'dir'); mkdir(outDir); end
nSamples = 200;

%% ========================== LOAD RETARDANCE STACK ==========================
retFiles = dir(fullfile(base_dir, retardanceGlob));
retFiles = retFiles(arrayfun(@(d) ~d.isdir && ~isempty(d.name) && ...
                                  d.name(1) ~= '.' && ...
                                  (isfield(d,'bytes') && d.bytes > 0), retFiles));
if isempty(retFiles)
    error('av_axis_analysis: no files matching %s in %s', retardanceGlob, base_dir);
end
nFramesRet = numel(retFiles);

%% ========================== PIV CELLS (from workspace) ==========================
assert(exist('u_original','var') && exist('v_original','var') && ...
       exist('x','var') && exist('y','var'), ...
    'av_axis_analysis: u_original, v_original, x, y must be in the workspace (run PIVlab first).');
Uc = u_original; Vc = v_original; Xc = x; Yc = y;
nFrames = min(nFramesRet, numel(Uc));
fprintf('av_axis_analysis: retardance=%d frames, PIV=%d frames, using=%d\n', ...
        nFramesRet, numel(Uc), nFrames);

%% ========================== CLICK AV ENDPOINTS ==========================
I0 = double(imread(fullfile(retFiles(1).folder, retFiles(1).name)));
fh = figure('Name','Click P1 then P2 (AV axis endpoints)');
imshow(I0, []); hold on;
% Optional: overlay the boundary polygon from a prior boundary_flows run.
if exist('resultsFile','var') && ~isempty(resultsFile) && exist(resultsFile,'file')==2
    Sres = load(resultsFile, 'visPolySeq', 'centroidXY', 'qcFlag');
    if isfield(Sres,'visPolySeq') && ~isempty(Sres.visPolySeq)
        fq = find(Sres.qcFlag, 1, 'first');
        if ~isempty(fq) && ~isempty(Sres.visPolySeq{fq})
            plot(Sres.visPolySeq{fq}(:,1), Sres.visPolySeq{fq}(:,2), ...
                 'c-', 'LineWidth', 1.2);
        end
    end
end
title('Click two points: P1 (animal pole) then P2 (vegetal pole)');
[xc, yc] = ginput(2);
assert(numel(xc) == 2, 'av_axis_analysis: need exactly 2 clicks.');
P1 = [xc(1), yc(1)];
P2 = [xc(2), yc(2)];
plot([P1(1) P2(1)], [P1(2) P2(2)], 'y-', 'LineWidth', 2);
plot(P1(1), P1(2), 'go', 'MarkerFaceColor','g', 'MarkerSize', 8);
plot(P2(1), P2(2), 'ro', 'MarkerFaceColor','r', 'MarkerSize', 8);
drawnow; pause(0.4); close(fh);

axis_vec       = P2 - P1;
axis_length_px = hypot(axis_vec(1), axis_vec(2));
d_AV           = axis_vec / axis_length_px;
d_perp         = [-d_AV(2), d_AV(1)];
s              = linspace(0, 1, nSamples)';
samplePts_px   = P1 + s * axis_vec;
samplePositions_um = (s * axis_length_px) / px_per_um;
fprintf('AV axis: P1=(%.1f, %.1f) -> P2=(%.1f, %.1f), length=%.1f px = %.2f um\n', ...
        P1(1),P1(2), P2(1),P2(2), axis_length_px, axis_length_px/px_per_um);

%% ========================== PER-FRAME SAMPLING ==========================
retardance_kymo = nan(nFrames, nSamples);
flow_along_kymo = nan(nFrames, nSamples);
flow_perp_kymo  = nan(nFrames, nSamples);

for fr = 1:nFrames
    % -- Retardance
    Iret = double(imread(fullfile(retFiles(fr).folder, retFiles(fr).name)));
    retardance_kymo(fr, :) = interp2(Iret, samplePts_px(:,1), samplePts_px(:,2), ...
                                     'linear', NaN);

    % -- PIV: scatteredInterpolant on the (possibly non-uniform) grid.
    % Cast to double — PIVlab can return single-precision arrays which
    % scatteredInterpolant refuses.
    X = double(Xc{fr}); Y = double(Yc{fr});
    U = double(Uc{fr}); V = double(Vc{fr});
    if isempty(X) || isempty(U), continue; end
    % Drop any NaN PIV samples before building the interpolant.
    ok = isfinite(X(:)) & isfinite(Y(:)) & isfinite(U(:)) & isfinite(V(:));
    if nnz(ok) < 4, continue; end
    Fu = scatteredInterpolant(X(ok), Y(ok), U(ok), 'linear', 'none');
    Fv = scatteredInterpolant(X(ok), Y(ok), V(ok), 'linear', 'none');
    u_s = Fu(samplePts_px(:,1), samplePts_px(:,2));
    v_s = Fv(samplePts_px(:,1), samplePts_px(:,2));

    % -- Convert to um/s
    switch lower(strtrim(pivVelUnit))
        case 'px_per_frame'
            u_s = (u_s / px_per_um) / dt_sec;
            v_s = (v_s / px_per_um) / dt_sec;
        case 'px_per_sec'
            u_s = u_s / px_per_um;
            v_s = v_s / px_per_um;
        otherwise
            error('av_axis_analysis: unknown pivVelUnit: %s', pivVelUnit);
    end

    flow_along_kymo(fr, :) = u_s * d_AV(1)   + v_s * d_AV(2);
    flow_perp_kymo(fr, :)  = u_s * d_perp(1) + v_s * d_perp(2);

    if mod(fr, 25) == 0
        fprintf('  processed frame %d/%d\n', fr, nFrames);
    end
end

time_min = (0:nFrames-1).' * (dt_sec/60);

%% ========================== SAVE ==========================
outMat = fullfile(outDir, 'av_axis_results.mat');
save(outMat, 'P1','P2','d_AV','d_perp','axis_length_px', ...
     'samplePts_px','samplePositions_um', ...
     'retardance_kymo','flow_along_kymo','flow_perp_kymo', ...
     'time_min','px_per_um','dt_sec','pivVelUnit','retardanceGlob','base_dir');
fprintf('Saved %s\n', outMat);

%% ========================== SUMMARY FIGURE ==========================
figS = figure('Position', [100 100 1500 500]);

subplot(1,3,1);
imagesc(samplePositions_um, time_min, retardance_kymo);
colormap(gca, 'hot'); cb = colorbar; ylabel(cb, 'Retardance (nm)');
xlabel('Position along AV (\mum, 0=P1, end=P2)'); ylabel('Time (min)');
title('Retardance along AV');

flim = max(abs([flow_along_kymo(:); flow_perp_kymo(:)]), [], 'omitnan');
if isempty(flim) || flim == 0 || ~isfinite(flim); flim = 1; end

subplot(1,3,2);
imagesc(samplePositions_um, time_min, flow_along_kymo); clim([-flim flim]);
colormap(gca, redblue_local(256)); cb = colorbar; ylabel(cb, 'v_{AV} (\mum/s, + toward P2)');
xlabel('Position along AV (\mum)'); ylabel('Time (min)');
title('Along-AV flow');

subplot(1,3,3);
imagesc(samplePositions_um, time_min, flow_perp_kymo); clim([-flim flim]);
colormap(gca, redblue_local(256)); cb = colorbar; ylabel(cb, 'v_{\perp} (\mum/s, + 90\circ CCW)');
xlabel('Position along AV (\mum)'); ylabel('Time (min)');
title('Perp-to-AV flow');

exportgraphics(figS, fullfile(outDir, 'av_axis_summary.png'), 'Resolution', 250);
fprintf('Saved %s\n', fullfile(outDir, 'av_axis_summary.png'));

%% ========================== PROFILES AT TIMEPOINTS ==========================
% Overlaid 1D profiles at nSnap equally-spaced times. Shows how the
% internal retardance and the internal flow components vary spatially
% along AV, at several moments through the movie.
nSnap = min(8, nFrames);
snapIdx = unique(round(linspace(1, nFrames, nSnap)));
cmapT = parula(numel(snapIdx));

figP = figure('Position', [100 100 1500 800]);

subplot(3,1,1); hold on; grid on;
for ii = 1:numel(snapIdx)
    plot(samplePositions_um, retardance_kymo(snapIdx(ii), :), ...
         'Color', cmapT(ii,:), 'LineWidth', 1.4);
end
ylabel('Retardance (nm)');
title('Internal retardance along AV at selected time points');

subplot(3,1,2); hold on; grid on;
for ii = 1:numel(snapIdx)
    plot(samplePositions_um, flow_along_kymo(snapIdx(ii), :), ...
         'Color', cmapT(ii,:), 'LineWidth', 1.4);
end
yline(0, 'k-');
ylabel('v_{AV} (\mum/s, + toward P2)');
title('Along-AV flow at selected time points');

subplot(3,1,3); hold on; grid on;
for ii = 1:numel(snapIdx)
    plot(samplePositions_um, flow_perp_kymo(snapIdx(ii), :), ...
         'Color', cmapT(ii,:), 'LineWidth', 1.4);
end
yline(0, 'k-');
ylabel('v_{\perp} (\mum/s, + 90\circ CCW)');
xlabel('Position along AV (\mum, 0=P1, end=P2)');
title('Perp-to-AV flow at selected time points');

% Shared time colorbar across the figure
cb = colorbar(gca, 'Position', [0.93 0.11 0.012 0.815]);
colormap(gca, cmapT);
clim(gca, [time_min(snapIdx(1)) time_min(snapIdx(end))]);
ylabel(cb, 'Time (min)');

exportgraphics(figP, fullfile(outDir, 'av_axis_profiles.png'), 'Resolution', 250);
fprintf('Saved %s\n', fullfile(outDir, 'av_axis_profiles.png'));

%% ========================== local helper ==========================
function cmap = redblue_local(n)
if nargin < 1; n = 256; end
half = floor(n/2);
low  = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
high = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
cmap = [low; high];
end

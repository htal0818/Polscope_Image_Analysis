function plot_wave_axis(resultsFile, outDir)
% PLOT_WAVE_AXIS  Visualize wave-axis results produced by wave_axis_pipeline.m.
%
% Produces four figures:
%   (a) kymograph_axis_overlay.png — Vtheta_kymo with phi_shape(t) and
%       theta_nuc(t) overlaid as lines.
%   (b) dphi_rose.png — rose plot of dPhi_shape for detected wave events.
%   (c) per_frame_axis.png — an example frame with boundary, centroid,
%       nucleus, and axis vector drawn.
%   (d) dphi_timeseries.png — time series of phi_shape, phi_vel, theta_nuc,
%       and dPhi with wave event markers.

if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(resultsFile), 'wave_axis');
end
if ~exist(outDir, 'dir'); mkdir(outDir); end

S = load(resultsFile);

nFrames = size(S.Vtheta_kymo, 1);
t_min   = S.time_min(:);

%% (a) Kymograph overlay
fh = figure('Visible','off','Position',[100 100 1100 650]);
imagesc(rad2deg(S.thetaCenters), t_min, S.Vtheta_kymo);
axis tight; colormap(parula);
cb = colorbar; ylabel(cb, 'v_\theta');
xlabel('Angle (deg)'); ylabel('Time (min)');
title('v_\theta(\theta,t) with wave axis (white) and nucleus direction (magenta)');
hold on;
% phi_shape(t) — unwrap to avoid jumps across 0/2pi
phi_plot = rad2deg(S.phi_shape);
nuc_plot = rad2deg(S.theta_nuc);
plot(phi_plot, t_min, 'w-',  'LineWidth', 1.5);
plot(nuc_plot, t_min, 'm--', 'LineWidth', 1.5);

% Mark detected events
if isfield(S, 'waveEvents') && ~isempty(S.waveEvents)
    fr_ev = [S.waveEvents.frame];
    fr_ev = fr_ev(fr_ev >= 1 & fr_ev <= nFrames);
    scatter(rad2deg(S.phi_shape(fr_ev)), t_min(fr_ev), 40, 'r', 'filled', ...
        'MarkerEdgeColor', 'w');
end
legend({'\phi_{shape}(t)', '\theta_{nuc}(t)', 'wave events'}, ...
       'Location', 'northeastoutside');
exportgraphics(fh, fullfile(outDir, 'kymograph_axis_overlay.png'), 'Resolution', 250);
close(fh);

%% (b) Rose plot of dPhi_shape for events
fh = figure('Visible','off','Position',[100 100 600 600]);
if isfield(S, 'waveEvents') && ~isempty(S.waveEvents)
    dphi = [S.waveEvents.dPhi_shape];
    dphi = dphi(isfinite(dphi));
    if ~isempty(dphi)
        polarhistogram(dphi, 18, 'FaceColor', [0.2 0.5 0.8]);
        title(sprintf('\\Delta\\phi (wave axis - nucleus direction), n=%d events\nCirc mean = %.1f\\circ,  R = %.3f,  Rayleigh p ~ %.2g', ...
            numel(dphi), rad2deg(S.circ_mean), S.circ_R, S.rayleigh_p));
    end
else
    text(0.5, 0.5, 'No events detected', 'HorizontalAlignment', 'center');
    axis off;
end
exportgraphics(fh, fullfile(outDir, 'dphi_rose.png'), 'Resolution', 250);
close(fh);

%% (c) Per-frame snapshot with boundary, centroid, nucleus, axis
% Pick the event with max amp_shape (most illustrative)
if isfield(S, 'waveEvents') && ~isempty(S.waveEvents)
    [~, idxBest] = max([S.waveEvents.amp_shape]);
    fr = S.waveEvents(idxBest).frame;
else
    fr = find(S.qcFlag, 1, 'first');
end

if ~isempty(fr) && fr <= nFrames
    fh = figure('Visible','off','Position',[100 100 700 700]);
    hold on; axis equal; axis ij;  % image coords: y down
    if isfield(S, 'visPolySeq') && ~isempty(S.visPolySeq{fr})
        poly = S.visPolySeq{fr};
        plot(poly(:,1), poly(:,2), 'k-', 'LineWidth', 1.5);
    end
    cxy = S.centroidXY(fr, :);
    nxy = S.nucleusXY(fr, :);
    plot(cxy(1), cxy(2), 'cs', 'MarkerFaceColor', 'c', 'MarkerSize', 10);
    text(cxy(1)+5, cxy(2), 'centroid');
    plot(nxy(1), nxy(2), 'm*', 'MarkerSize', 14, 'LineWidth', 2);
    text(nxy(1)+5, nxy(2), 'nucleus');

    % Axis vector from centroid
    phi = S.phi_shape(fr);
    Lax = hypot(nxy(1)-cxy(1), nxy(2)-cxy(2)) + 30;
    quiver(cxy(1), cxy(2), Lax*cos(phi), Lax*sin(phi), 0, ...
        'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    legend({'boundary','centroid','nucleus','\phi_{shape}'}, 'Location','best');
    title(sprintf('Frame %d: wave axis (red) vs centroid->nucleus', fr));
    exportgraphics(fh, fullfile(outDir, sprintf('per_frame_axis_fr%04d.png', fr)), ...
                   'Resolution', 250);
    close(fh);
end

%% (d) Time series
fh = figure('Visible','off','Position',[100 100 1100 700]);

subplot(3,1,1);
plot(t_min, rad2deg(S.phi_shape), 'b-', 'LineWidth', 1.2); hold on;
plot(t_min, rad2deg(S.phi_vel),   'g-', 'LineWidth', 1.0);
plot(t_min, rad2deg(S.theta_nuc), 'm--','LineWidth', 1.2);
ylabel('Angle (deg)'); grid on;
legend({'\phi_{shape}','\phi_{vel}','\theta_{nuc}'}, 'Location', 'best');
title('Wave axes vs. nucleus direction over time');

subplot(3,1,2);
plot(t_min, rad2deg(S.dPhi_shape), 'b-', 'LineWidth', 1.2); hold on;
plot(t_min, rad2deg(S.dPhi_vel),   'g-', 'LineWidth', 1.0);
yline(0, 'k--');
ylabel('\Delta\phi (deg)'); grid on;
legend({'shape - nuc','vel - nuc'}, 'Location', 'best');

subplot(3,1,3);
plot(t_min, S.amp_shape, 'b-', 'LineWidth', 1.2); hold on;
plot(t_min, S.amp_vel,   'g-', 'LineWidth', 1.0);
if isfield(S, 'waveEvents') && ~isempty(S.waveEvents)
    fr_ev = [S.waveEvents.frame];
    fr_ev = fr_ev(fr_ev >= 1 & fr_ev <= nFrames);
    scatter(t_min(fr_ev), S.amp_shape(fr_ev), 30, 'r', 'filled');
end
xlabel('Time (min)'); ylabel('Amplitude'); grid on;
legend({'|c_{k*}^{shape}|','|c_1^{v}|','events'}, 'Location', 'best');

exportgraphics(fh, fullfile(outDir, 'dphi_timeseries.png'), 'Resolution', 250);
close(fh);

end

function plot_wave_axis(resultsFile, outDir)
% PLOT_WAVE_AXIS  Visualize wave-axis results produced by wave_axis_pipeline.m.
%
% Gracefully handles runs where the shape-axis (phi_shape / amp_shape) is
% all NaN because visPolySeq wasn't in the source .mat — falls back to
% showing the velocity-axis (phi_vel / amp_vel) as the primary metric.
%
% Produces four figures in outDir:
%   (a) kymograph_axis_overlay.png
%   (b) dphi_rose.png
%   (c) per_frame_axis_fr####.png
%   (d) dphi_timeseries.png

if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(resultsFile), 'wave_axis');
end
if ~exist(outDir, 'dir'); mkdir(outDir); end

S = load(resultsFile);

nFrames = size(S.Vtheta_kymo, 1);
t_min   = S.time_min(:);

% --- Decide which axis is "primary": shape if we have any finite values,
% otherwise fall back to vel. This mirrors wave_axis_pipeline.m.
if isfield(S, 'phi_shape') && any(isfinite(S.phi_shape))
    phi_primary   = S.phi_shape;
    amp_primary   = S.amp_shape;
    dPhi_primary  = S.dPhi_shape;
    primaryLabel  = '\phi_{shape}';
    ampLabel      = '|c_{k*}^{shape}|';
    dphiField     = 'dPhi_shape';
else
    phi_primary   = S.phi_vel;
    amp_primary   = S.amp_vel;
    dPhi_primary  = S.dPhi_vel;
    primaryLabel  = '\phi_{vel}';
    ampLabel      = '|c_1^{v}|';
    dphiField     = 'dPhi_vel';
end

%% (a) Kymograph overlay
fh = figure('Visible','off','Position',[100 100 1100 650]);
imagesc(rad2deg(S.thetaCenters), t_min, S.Vtheta_kymo);
axis tight; colormap(parula);
cb = colorbar; ylabel(cb, 'v_\theta');
xlabel('Angle (deg)'); ylabel('Time (min)');
title(sprintf('v_\\theta(\\theta,t) with %s (white) and \\theta_{nuc} (magenta)', primaryLabel));
hold on;
plot(rad2deg(phi_primary), t_min, 'w-',  'LineWidth', 1.5);
plot(rad2deg(S.theta_nuc), t_min, 'm--', 'LineWidth', 1.5);

if isfield(S, 'waveEvents') && ~isempty(S.waveEvents)
    fr_ev = [S.waveEvents.frame];
    fr_ev = fr_ev(fr_ev >= 1 & fr_ev <= nFrames);
    scatter(rad2deg(phi_primary(fr_ev)), t_min(fr_ev), 40, 'r', 'filled', ...
        'MarkerEdgeColor', 'w');
end
legend({primaryLabel, '\theta_{nuc}(t)', 'wave events'}, ...
       'Location', 'northeastoutside');
exportgraphics(fh, fullfile(outDir, 'kymograph_axis_overlay.png'), 'Resolution', 250);
close(fh);

%% (b) Rose plot of dPhi for events
fh = figure('Visible','off','Position',[100 100 600 600]);
if isfield(S, 'waveEvents') && ~isempty(S.waveEvents) && isfield(S.waveEvents, dphiField)
    dphi = arrayfun(@(e) e.(dphiField), S.waveEvents);
    dphi = dphi(isfinite(dphi));
    if ~isempty(dphi)
        polarhistogram(dphi, 18, 'FaceColor', [0.2 0.5 0.8]);
        title(sprintf('\\Delta\\phi (%s - \\theta_{nuc}), n=%d events\nCirc mean = %.1f\\circ,  R = %.3f,  Rayleigh p ~ %.2g', ...
            primaryLabel, numel(dphi), rad2deg(S.circ_mean), S.circ_R, S.rayleigh_p));
    end
else
    text(0.5, 0.5, 'No events detected', 'HorizontalAlignment', 'center');
    axis off;
end
exportgraphics(fh, fullfile(outDir, 'dphi_rose.png'), 'Resolution', 250);
close(fh);

%% (c) Per-frame snapshot with boundary (if any), centroid, nucleus, axis
if isfield(S, 'waveEvents') && ~isempty(S.waveEvents)
    % Pick the strongest-amplitude event, using whichever amplitude the
    % events stored. waveEvents added in the new pipeline has .amp;
    % older files may have .amp_shape.
    if isfield(S.waveEvents, 'amp')
        ampArr = [S.waveEvents.amp];
    else
        ampArr = [S.waveEvents.amp_shape];
    end
    [~, idxBest] = max(ampArr);
    fr = S.waveEvents(idxBest).frame;
else
    fr = find(S.qcFlag, 1, 'first');
end

if ~isempty(fr) && fr <= nFrames
    fh = figure('Visible','off','Position',[100 100 700 700]);
    hold on; axis equal; axis ij;
    if isfield(S, 'visPolySeq') && ~isempty(S.visPolySeq) && ...
            ~isempty(S.visPolySeq{fr})
        poly = S.visPolySeq{fr};
        plot(poly(:,1), poly(:,2), 'k-', 'LineWidth', 1.5);
        havePoly = true;
    else
        havePoly = false;
    end
    cxy = S.centroidXY(fr, :);
    nxy = S.nucleusXY(fr, :);
    plot(cxy(1), cxy(2), 'cs', 'MarkerFaceColor', 'c', 'MarkerSize', 10);
    text(cxy(1)+5, cxy(2), 'centroid');
    plot(nxy(1), nxy(2), 'm*', 'MarkerSize', 14, 'LineWidth', 2);
    text(nxy(1)+5, nxy(2), 'nucleus');

    phi = phi_primary(fr);
    if isfinite(phi)
        Lax = hypot(nxy(1)-cxy(1), nxy(2)-cxy(2)) + 30;
        quiver(cxy(1), cxy(2), Lax*cos(phi), Lax*sin(phi), 0, ...
            'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    end
    if havePoly
        legend({'boundary','centroid','nucleus',primaryLabel}, 'Location','best');
    else
        legend({'centroid','nucleus',primaryLabel}, 'Location','best');
    end
    title(sprintf('Frame %d: wave axis (red) vs centroid->nucleus', fr));
    exportgraphics(fh, fullfile(outDir, sprintf('per_frame_axis_fr%04d.png', fr)), ...
                   'Resolution', 250);
    close(fh);
end

%% (d) Time series
fh = figure('Visible','off','Position',[100 100 1100 700]);

subplot(3,1,1);
if any(isfinite(S.phi_shape))
    plot(t_min, rad2deg(S.phi_shape), 'b-', 'LineWidth', 1.2); hold on;
end
plot(t_min, rad2deg(S.phi_vel),   'g-', 'LineWidth', 1.0); hold on;
plot(t_min, rad2deg(S.theta_nuc), 'm--','LineWidth', 1.2);
ylabel('Angle (deg)'); grid on;
if any(isfinite(S.phi_shape))
    legend({'\phi_{shape}','\phi_{vel}','\theta_{nuc}'}, 'Location', 'best');
else
    legend({'\phi_{vel}','\theta_{nuc}'}, 'Location', 'best');
end
title('Wave axes vs. nucleus direction over time');

subplot(3,1,2);
if any(isfinite(S.dPhi_shape))
    plot(t_min, rad2deg(S.dPhi_shape), 'b-', 'LineWidth', 1.2); hold on;
end
plot(t_min, rad2deg(S.dPhi_vel), 'g-', 'LineWidth', 1.0); hold on;
yline(0, 'k--');
ylabel('\Delta\phi (deg)'); grid on;
if any(isfinite(S.dPhi_shape))
    legend({'shape - nuc','vel - nuc'}, 'Location', 'best');
else
    legend({'vel - nuc'}, 'Location', 'best');
end

subplot(3,1,3);
if any(isfinite(S.amp_shape))
    plot(t_min, S.amp_shape, 'b-', 'LineWidth', 1.2); hold on;
end
plot(t_min, S.amp_vel, 'g-', 'LineWidth', 1.0); hold on;
if isfield(S, 'waveEvents') && ~isempty(S.waveEvents)
    fr_ev = [S.waveEvents.frame];
    fr_ev = fr_ev(fr_ev >= 1 & fr_ev <= nFrames);
    scatter(t_min(fr_ev), amp_primary(fr_ev), 30, 'r', 'filled');
end
xlabel('Time (min)'); ylabel('Amplitude'); grid on;
if any(isfinite(S.amp_shape))
    legend({'|c_{k*}^{shape}|','|c_1^{v}|','events'}, 'Location', 'best');
else
    legend({ampLabel,'events'}, 'Location', 'best');
end

exportgraphics(fh, fullfile(outDir, 'dphi_timeseries.png'), 'Resolution', 250);
close(fh);

end

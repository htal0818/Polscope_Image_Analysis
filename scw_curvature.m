function C = scw_curvature(R, varargin)
%SCW_CURVATURE  Curvature change and SCW strength from SCW_CORTEX output.
%
%   C = SCW_CURVATURE(R) quantifies the surface contraction wave (SCW) from
%   the per-frame contours in R (output of SCW_CORTEX), following the
%   curvature schema of Bischof et al., Nat Commun 8:849 (2017), adapted to
%   LC-PolScope retardance in place of fluorescence:
%
%     1. SMOOTH   The segmented outline of each frame is expressed in polar
%                 coordinates r(theta) about the frame centroid and smoothed
%                 by piecewise polynomial fits: an order-'PolyOrder' local
%                 least-squares fit over a sliding window of 'SmoothUm'
%                 micrometres of arc (circular Savitzky-Golay). Where the
%                 paper segments with Chan-Vese, SCW_CORTEX's ridge-snapped
%                 contour is used instead - it sits on the retardance
%                 maximum of the cortex with 0.086 px rms error.
%     2. CURVE    The first principal curvature is computed for small
%                 segments of the outline, each about 'SegUm' (2 um) of
%                 arc, from derivative stencils of the local fit (the
%                 least-squares generalisation of finite differences):
%                    kappa = (r^2 + 2 r'^2 - r r'') / (r^2 + r'^2)^(3/2)
%                 kappa = 1/r for a circle; negative in concave dents.
%     3. RELATIVE To remove uneven starting shapes, relative curvature is
%                 kappa minus the curvature of the reference frame
%                 ('RefFrames', default the first QC-ok frame), segment by
%                 segment.
%     4. STRENGTH SCW strength = variance of the radii of curvature during
%                 'ScwWindow' minus the variance during 'BgWindow', an
%                 equal time window in metaphase. Radii are capped at
%                 'RhoCapUm' so near-flat segments cannot dominate.
%     5. RINGS    Cortical and subcortical retardance (R.cortNm, R.subNm -
%                 the rings SCW_CORTEX builds around the segmented outline)
%                 are averaged in the SAME theta segments used for the
%                 curvature, so curvature and cortex signal can be compared
%                 segment for segment.
%
%   The segment grid is fixed from the first usable frame (M segments of
%   ~SegUm arc at theta = 0, 360/M, 2*360/M, ... degrees, theta measured
%   from the +x axis toward +y/image-down about the per-frame centroid), so
%   rows of every output matrix align across time and kymographs are
%   registered. The outline is assumed star-shaped about its centroid
%   (r(theta) single-valued), which holds for oocytes.
%
%   Options (name-value)
%     'SegUm'      arc length per curvature segment, um     default 2
%     'SmoothUm'   arc window of the piecewise poly fit, um default 16
%     'PolyOrder'  order of the local polynomial (>= 2)     default 3
%     'RefFrames'  frame number(s) for relative curvature   default first ok
%     'ScwWindow'  [t0 t1] min containing the SCW           default [] (skip)
%     'BgWindow'   [t0 t1] min in metaphase, same duration  default [] (auto)
%     'RhoCapUm'   cap on |radius of curvature|, um         default 250
%     'OutDir'     write CSVs here                          default '' (skip)
%     'Plot'       summary figure                           default true
%
%   Output C (nF frames x M segments unless noted):
%     segThetaDeg  1 x M   segment centre angles (deg)
%     segUm        scalar  nominal arc length per segment (um)
%     rUm          smoothed radius per segment (um)
%     kappa        curvature (1/um)
%     kappaRel     curvature minus reference frame (1/um)
%     rho          radius of curvature, signed, capped (um)
%     cortNm       mean cortical retardance per segment (nm)
%     subNm        mean subcortical retardance per segment (nm)
%     rhoVarTime   nF x 1  var of rho across segments per frame (um^2)
%     scw          struct: window, bgWindow, varScw, varBg, strengthUm2,
%                  strengthSegUm2 (1 x M, per-segment temporal variance
%                  difference), nScw, nBg
%     frames, timeMin, ok, refFrames, smoothUmUsed
%
%   Typical use:
%     R = scw_cortex(src);
%     C = scw_curvature(R);                        % look at rhoVarTime
%     C = scw_curvature(R, 'ScwWindow',[30 55], 'BgWindow',[0 25], ...
%                          'OutDir', fullfile(src,'scw_output'));

p = inputParser;
p.addParameter('SegUm',      2);
p.addParameter('SmoothUm',  16);
p.addParameter('PolyOrder',  3);
p.addParameter('RefFrames', []);
p.addParameter('ScwWindow', []);
p.addParameter('BgWindow',  []);
p.addParameter('RhoCapUm', 250);
p.addParameter('OutDir',   '');
p.addParameter('Plot',     true);
p.parse(varargin{:});
o = p.Results;
if o.PolyOrder < 2
    error('scw:polyOrder', 'PolyOrder must be >= 2 to give a second derivative.');
end

nF = numel(R.frames);
haveSnake = ~cellfun(@isempty, R.snake);
k0 = find(haveSnake & R.ok, 1);
if isempty(k0), k0 = find(haveSnake, 1); end
if isempty(k0), error('scw:noContours', 'R contains no contours.'); end

% ---- fixed segment grid from the reference geometry
rbar0 = mean(hypot(R.snake{k0}(:,1) - R.centre(k0,1), ...
                   R.snake{k0}(:,2) - R.centre(k0,2))) / R.pxPerUm;
M   = max(16, round(2*pi*rbar0 / o.SegUm));    % segments of ~SegUm arc
Nf  = 8*M;                                     % fine grid, 8 per segment
thf = (0:Nf-1)' * (2*pi/Nf);
dthf = 2*pi/Nf;
ths = thf(1:8:end);                            % segment centres = subsample

% ---- piecewise-polynomial (Savitzky-Golay) stencils on the fine grid.
% The window is SmoothUm of arc; the derivative rows of the fit are the
% least-squares finite-difference stencils used for the curvature.
arcPt = 2*pi*rbar0 / Nf;                       % um of arc per fine point
win = round(o.SmoothUm / arcPt);
if mod(win,2) == 0, win = win + 1; end
win = max(win, 2*floor((o.PolyOrder+2)/2) + 1);
half = (win-1)/2;
xw = (-half:half)';
V  = xw .^ (0:o.PolyOrder);
P  = pinv(V);
k0v = P(1,:)';  k1v = P(2,:)';  k2v = 2*P(3,:)';

hasCort = isfield(R, 'cortNm');
if ~hasCort
    warning('scw:noRings', ['R has no cortNm/subNm (older scw_cortex). ' ...
        'Cortical signal falls back to R.kymo (nm*px integral); ' ...
        'subcortical signal will be NaN.']);
end

C = struct();
C.segThetaDeg = ths' * 180/pi;
C.segUm   = 2*pi*rbar0 / M;
C.smoothUmUsed = win * arcPt;
C.frames  = R.frames;  C.timeMin = R.timeMin;  C.ok = R.ok;
[C.rUm, C.kappa, C.cortNm, C.subNm] = deal(nan(nF, M));

for k = 1:nF
    if ~haveSnake(k), continue; end
    sk  = R.snake{k};  cen = R.centre(k,:);
    thP = mod(atan2(sk(:,1)-cen(1), sk(:,2)-cen(2)), 2*pi);  % per point
    rP  = hypot(sk(:,1)-cen(1), sk(:,2)-cen(2)) / R.pxPerUm;

    % r(theta) on the fine grid, circular
    [th, is] = sort(thP);  r = rP(is);
    keep = [true; diff(th) > 1e-9];  th = th(keep);  r = r(keep);
    if numel(th) < 8, continue; end
    rf = interp1([th-2*pi; th; th+2*pi], [r; r; r], thf, 'linear');

    r0 = circCorr(rf, k0v);
    r1 = circCorr(rf, k1v) / dthf;
    r2 = circCorr(rf, k2v) / dthf^2;
    kf = (r0.^2 + 2*r1.^2 - r0.*r2) ./ (r0.^2 + r1.^2).^1.5;

    C.rUm(k,:)   = r0(1:8:end)';
    C.kappa(k,:) = kf(1:8:end)';

    % ring signal, averaged in the same theta segments
    bi = mod(round(thP * M/(2*pi)), M) + 1;
    if hasCort
        C.cortNm(k,:) = accumarray(bi, R.cortNm(k,:)', [M 1], ...
                                   @(v) mean(v,'omitnan'), NaN)';
        C.subNm(k,:)  = accumarray(bi, R.subNm(k,:)',  [M 1], ...
                                   @(v) mean(v,'omitnan'), NaN)';
    else
        C.cortNm(k,:) = accumarray(bi, R.kymo(k,:)',   [M 1], ...
                                   @(v) mean(v,'omitnan'), NaN)';
    end
end

% ---- relative curvature: subtract the reference frame's curvature
if isempty(o.RefFrames)
    refIdx = k0;
else
    refIdx = find(ismember(R.frames, o.RefFrames));
    if isempty(refIdx)
        error('scw:refFrames', 'RefFrames not found among processed frames.');
    end
end
C.refFrames = R.frames(refIdx);
C.kappaRel  = C.kappa - mean(C.kappa(refIdx,:), 1, 'omitnan');

% ---- radii of curvature, capped, and their spread
C.rho = 1 ./ C.kappa;
big = abs(C.rho) > o.RhoCapUm;
C.rho(big) = sign(C.rho(big)) * o.RhoCapUm;
C.rhoVarTime = var(C.rho, 0, 2, 'omitnan');

% ---- SCW strength: variance of radii of curvature, background-subtracted
C.scw = struct('window', o.ScwWindow, 'bgWindow', o.BgWindow, ...
               'varScw', NaN, 'varBg', NaN, 'strengthUm2', NaN, ...
               'strengthSegUm2', nan(1, M), 'nScw', 0, 'nBg', 0);
if ~isempty(o.ScwWindow)
    bg = o.BgWindow;
    if isempty(bg)
        bg = [R.timeMin(1), R.timeMin(1) + diff(o.ScwWindow)];
        warning('scw:bgWindow', ['No BgWindow given; using the first %.1f ' ...
            'min of the recording as the metaphase background.'], diff(bg));
    end
    if abs(diff(bg) - diff(o.ScwWindow)) > 0.1 * diff(o.ScwWindow)
        warning('scw:bgWindow', ['BgWindow (%.1f min) and ScwWindow ' ...
            '(%.1f min) durations differ; the variance comparison is ' ...
            'biased by the unequal sample counts.'], diff(bg), diff(o.ScwWindow));
    end
    inS = R.timeMin >= o.ScwWindow(1) & R.timeMin <= o.ScwWindow(2) & R.ok;
    inB = R.timeMin >= bg(1)          & R.timeMin <= bg(2)          & R.ok;
    C.scw.bgWindow = bg;
    C.scw.nScw = nnz(inS);  C.scw.nBg = nnz(inB);
    if nnz(inS) < 2 || nnz(inB) < 2
        warning('scw:window', 'Fewer than 2 QC-ok frames in a window; strength is NaN.');
    else
        vS = C.rho(inS,:);  vB = C.rho(inB,:);
        C.scw.varScw = var(vS(:), 0, 'omitnan');
        C.scw.varBg  = var(vB(:), 0, 'omitnan');
        C.scw.strengthUm2 = C.scw.varScw - C.scw.varBg;
        C.scw.strengthSegUm2 = var(vS, 0, 1, 'omitnan') - var(vB, 0, 1, 'omitnan');
    end
end

fprintf(['curvature: %d segments of %.2f um arc, poly order %d over ' ...
         '%.1f um, reference frame(s) %s\n'], M, C.segUm, o.PolyOrder, ...
        C.smoothUmUsed, mat2str(C.refFrames(:)'));
if ~isempty(o.ScwWindow) && isfinite(C.scw.strengthUm2)
    fprintf(['SCW strength: var(rho) %.1f um^2 in [%g %g] min minus %.1f ' ...
             'um^2 background in [%g %g] min = %.1f um^2 ' ...
             '(%d and %d frames)\n'], C.scw.varScw, o.ScwWindow(1), ...
            o.ScwWindow(2), C.scw.varBg, C.scw.bgWindow(1), ...
            C.scw.bgWindow(2), C.scw.strengthUm2, C.scw.nScw, C.scw.nBg);
end

if ~isempty(o.OutDir), exportCsv(C, R, o); end
if o.Plot, plotCurvature(C, R, o); end
end


% ===================================================================== core

function out = circCorr(v, k)
%CIRCCORR  Circular correlation of column v with stencil k(-h..h).
h  = (numel(k)-1)/2;
vp = [v(end-h+1:end); v; v(1:h)];
out = conv(vp, flipud(k(:)), 'valid');
end


% ================================================================== outputs

function exportCsv(C, R, o)
if ~exist(o.OutDir, 'dir'), mkdir(o.OutDir); end
nF = numel(C.frames);  M = numel(C.segThetaDeg);

T = table(repelem(C.frames(:), M), repelem(C.timeMin(:), M), ...
          repmat((1:M)', nF, 1), repmat(C.segThetaDeg(:), nF, 1), ...
          reshape(C.rUm', [], 1), reshape(C.kappa', [], 1), ...
          reshape(C.kappaRel', [], 1), reshape(C.rho', [], 1), ...
          reshape(C.cortNm', [], 1), reshape(C.subNm', [], 1), ...
          repelem(C.ok(:), M), ...
    'VariableNames', {'frame','time_min','segment','theta_deg','r_um', ...
        'kappa_per_um','kappa_rel_per_um','rho_um','cortical_nm', ...
        'subcortical_nm','qc_ok'});
writetable(T, fullfile(o.OutDir, 'scw_curvature_segments.csv'));

F = table(C.frames(:), C.timeMin(:), C.rhoVarTime(:), ...
          mean(C.kappa, 2, 'omitnan'), mean(C.cortNm, 2, 'omitnan'), ...
          mean(C.subNm, 2, 'omitnan'), C.ok(:), ...
    'VariableNames', {'frame','time_min','rho_var_um2','kappa_mean_per_um', ...
        'cortical_mean_nm','subcortical_mean_nm','qc_ok'});
writetable(F, fullfile(o.OutDir, 'scw_curvature_frames.csv'));

if isfinite(C.scw.strengthUm2)
    W = table(C.scw.window(1), C.scw.window(2), C.scw.bgWindow(1), ...
              C.scw.bgWindow(2), C.scw.varScw, C.scw.varBg, ...
              C.scw.strengthUm2, C.scw.nScw, C.scw.nBg, ...
        'VariableNames', {'scw_t0_min','scw_t1_min','bg_t0_min','bg_t1_min', ...
            'var_scw_um2','var_bg_um2','strength_um2','n_frames_scw', ...
            'n_frames_bg'});
    writetable(W, fullfile(o.OutDir, 'scw_strength.csv'));
end
fprintf('wrote CSVs to %s\n', o.OutDir);
end


function plotCurvature(C, R, o)
t  = C.timeMin;  th = C.segThetaDeg;
um = 1/R.pxPerUm;

f = findobj('Type','figure','Tag','scw_curv');
if isempty(f)
    f = figure('Color','k','Position',[60 60 1500 820],'Tag','scw_curv');
else, figure(f); clf(f);
end
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% 1 -- reference outline coloured by curvature
nexttile;
kref = find(ismember(R.frames, C.refFrames), 1);
cen  = R.centre(kref,:);
xs = cen(2)*um + C.rUm(kref,:) .* cos(th*pi/180);
ys = cen(1)*um + C.rUm(kref,:) .* sin(th*pi/180);
scatter(xs, ys, 18, C.kappa(kref,:), 'filled');
axis image; set(gca,'YDir','reverse','Color','k','XColor','w','YColor','w');
colormap(gca, parula); cb = colorbar('Color','w'); cb.Label.String = '1/um';
xlabel('x (um)'); ylabel('y (um)');
title(sprintf('reference frame %d curvature', C.refFrames(1)), 'Color','w');

% 2 -- relative curvature kymograph
nexttile([1 2]);
K = C.kappaRel;
imagesc(th, t, K, 'AlphaData', ~isnan(K)); set(gca,'Color','k');
colormap(gca, divMap(256));
m = prctile(abs(K(~isnan(K))), 98);
if isfinite(m) && m > 0, clim([-m m]); end
set(gca,'XColor','w','YColor','w'); colorbar('Color','w');
xlabel('theta (deg)'); ylabel('time (min)');
title('relative curvature \Delta\kappa (1/um): red = locally contracted', 'Color','w');

% 3 -- spread of radii of curvature over time, with the analysis windows
nexttile;
plot(t, C.rhoVarTime, 'c', 'LineWidth', 1.4); hold on;
plot(t(~C.ok), C.rhoVarTime(~C.ok), 'r.', 'MarkerSize', 10);
yl = ylim;
shade = @(w, col) patch([w(1) w(2) w(2) w(1)], [yl(1) yl(1) yl(2) yl(2)], ...
    col, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
if ~isempty(C.scw.window),   shade(C.scw.window,   [1 .3 .3]); end
if ~isempty(C.scw.bgWindow), shade(C.scw.bgWindow, [.3 .5 1]); end
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('time (min)'); ylabel('var(\rho) across segments (um^2)');
if isfinite(C.scw.strengthUm2)
    title(sprintf('SCW strength = %.1f um^2', C.scw.strengthUm2), 'Color','w');
else
    title('var of radii of curvature (pick ScwWindow here)', 'Color','w');
end

% 4 -- cortical retardance in the same segments
nexttile([1 2]);
K = C.cortNm;
imagesc(th, t, K, 'AlphaData', ~isnan(K)); set(gca,'Color','k');
colormap(gca, hot);
cl = prctile(K(~isnan(K)), [2 98]);
if numel(cl) == 2 && cl(2) > cl(1), clim(cl); end
set(gca,'XColor','w','YColor','w'); colorbar('Color','w');
xlabel('theta (deg)'); ylabel('time (min)');
title('cortical retardance in curvature segments (nm)', 'Color','w');

% 5 -- subcortical retardance, or per-segment strength if computed
nexttile;
if any(isfinite(C.subNm(:)))
    K = C.subNm;
    imagesc(th, t, K, 'AlphaData', ~isnan(K)); set(gca,'Color','k');
    colormap(gca, hot);
    cl = prctile(K(~isnan(K)), [2 98]);
    if numel(cl) == 2 && cl(2) > cl(1), clim(cl); end
    set(gca,'XColor','w','YColor','w'); colorbar('Color','w');
    xlabel('theta (deg)'); ylabel('time (min)');
    title('subcortical retardance (nm)', 'Color','w');
elseif any(isfinite(C.scw.strengthSegUm2))
    plot(th, C.scw.strengthSegUm2, 'c', 'LineWidth', 1.4);
    set(gca,'Color','k','XColor','w','YColor','w');
    xlabel('theta (deg)'); ylabel('strength (um^2)');
    title('per-segment SCW strength', 'Color','w');
else
    axis off;
end
end


function cmap = divMap(n)
%DIVMAP  Blue-white-red diverging colormap (negative-zero-positive).
half = floor(n/2);
up = linspace(0, 1, half)';
cmap = [ [up, up, ones(half,1)]; ones(n-2*half, 3); ...
         [ones(half,1), flipud(up), flipud(up)] ];
end

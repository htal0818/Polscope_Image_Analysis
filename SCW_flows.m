close all;
%% ========================== USER INPUTS ==================================
% --- State-1 input: set ONE of these forms ---
%  A) Single file:
% state1_input = '/ABS/PATH/to/state1_frame.tif';

%  B) Folder of images:
% state1_input = '/ABS/PATH/to/state1_folder/';



% test this on frames where SCW is actively occuring

% find a maximum speed the flows hit

% profile on

state1_input = '/Users/hridaytalreja/Desktop/June_2026_data/2026_06_17_F11/SMS_2026_0617_1256_1/Pos0/*State1*';

% --- Output root folder ---
outDir = '/Users/hridaytalreja/Desktop/June_2026_data/2026_06_17_F11/SMS_2026_0617_1256_1/Pos0/PIV_Jun18/';
mkdir(outDir);
OvlerayDir = '/Users/hridaytalreja/Desktop/June_2026_data/2026_06_17_F11/SMS_2026_0617_1256_1/Pos0/PIV_Jun18/Overlays/';
mkdir(OvlerayDir);


% % % %


base_dir = '/Users/hridaytalreja/Desktop/June_2026_data/2026_06_17_F11/SMS_2026_0617_1256_1/Pos0';


s1 = strcat(base_dir,'/*State1*.tif');
s2 = strcat(base_dir,'/*State2*.tif');
s3 = strcat(base_dir,'/*State3*.tif');
s4 = strcat(base_dir,'/*State4*.tif');

d1 = dir(s1);
d2 = dir(s2);
d3 = dir(s3);
d4 = dir(s4);


% % %


px_per_um = 6.25/2;       % 6.25 px per um (40x objective) or 3.12 for 20x
dt        = 20;         % 30 s per frame

% --- Overlay / quiver appearance (single source of truth) ---
quiverColor   = [1 0 0];   % color of BOTH the flow arrows and the scale bar
arrowScale    = 10;        % velocity -> length: plot-pixels per (pixel/frame) of displacement
                           %   the ONE scale used for the arrows AND the bar
                           %   (increase if arrows look small; the bar rescales with it)
arrowStride   = 3;         % plot every Nth vector in x & y (higher = fewer arrows)
arrowHeadSize = 1.9;         % flow-arrow head size (MaxHeadSize)

% --- Velocity scale bar (reference arrow, upper-right corner) ---
scaleBar_refVel_umin = 5; % speed the reference arrow represents (um/min) = the label
scaleBar_backdrop    = true;

% --- Save overlays ---
makeVideo = true;

%% ================== PIVlab arrays from workspace =========================

Uc = u_original; Vc = v_original; Xc = x; Yc = y;

nFrames = numel(Uc);
t_sec = (0:nFrames-1) * dt;  %  conversion to s
t_min = t_sec/60; %conversion to min

% grid spacing (pixels) for divergence after unit conversion
dx_px = median(diff(unique(Xc{1}(:))));
dy_px = median(diff(unique(Yc{1}(:))));


%% =================== Segmentation params (contour_retardance pipeline) ===
segParams.sigmaBlur            = 1;       % Gaussian blur sigma (px)
segParams.closeRadius          = 25;      % morphological close disk radius (px) — matches contour_retardance.m
segParams.minArea              = 5000;    % min object area (px^2)
segParams.thresholdMode        = 'otsu';  % 'adaptive', 'edge', or 'gradient'
segParams.segFromMask          = true;    % true for 4-state (oocyte is dark)
segParams.adaptSensitivity     = 1;
segParams.adaptNeighborhood    = 201;
segParams.edgeMethod           = 'Sobel';
segParams.edgeDilateRadius     = 2;
segParams.gradientPercentile   = 80;
segParams.useCaching           = false;
segParams.cacheIntensityThreshold = 0.02;
segParams.cacheForceRecalcEveryN  = 10;

% --- Boundary smoothing / curvature diagnostics ---
% AC is useful when the thresholded edge is close but jagged; turn off if it
% drifts onto internal texture for a given oocyte/movie.
segParams.useActiveContour       = false;
segParams.activeContourIterations = 100;
segParams.activeContourMethod    = 'Chan-Vese';
segParams.smoothMask             = false;   % match contour.m: keep the raw Otsu mask (no rebuild from smoothed boundary)
segParams.smoothBoundary         = true;    % smooth only the RETURNED polygon (curvature/geom), not the mask
segParams.nBoundaryPts           = 720;
segParams.boundarySmoothFrac     = 0.05;
segParams.forceOuterEnvelope     = false;   % match contour.m: no outer-envelope rebuild; cleanup is imclose+imfill
segParams.outerEnvelopePercentile = 95;

segCache = [];
maskErodePx = 1;   % erode mask by N px before PIV masking (exclude edge vectors)

%% =================== Containers =========================================
BW_seq     = cell(nFrames,1);
polySeq    = cell(nFrames,1);
geomSeq    = cell(nFrames,1);
INgridSeq  = cell(nFrames,1);
U_masked   = cell(nFrames,1);
V_masked   = cell(nFrames,1);

nBoundaryPts = segParams.nBoundaryPts;
xc_px = nan(nFrames,1);
yc_px = nan(nFrames,1);
R_fit_px = nan(nFrames,1);
curvatureByArc_umInv = nan(nFrames, nBoundaryPts);
radiusCurvatureByArc_um = nan(nFrames, nBoundaryPts);
thetaByArc_deg = nan(nFrames, nBoundaryPts);
meanRadiusCurvature_um = nan(nFrames,1);
medianRadiusCurvature_um = nan(nFrames,1);

%% =================== MAIN LOOP: build mask (State-1) + map to PIV =======
for t = 1:nFrames


    maskFile = fullfile(outDir, sprintf('mask_%04d.png', t));
    overlayFile = fullfile(outDir, sprintf('overlay_%04d.png', t));
    plotFile = fullfile(outDir, 'velocity_vs_time.png');
    csvFile  = fullfile(outDir, 'flow_metrics.csv');

    a1  = im2double(imread(fullfile(d1(t).folder, d1(t).name)));
    a2 = im2double(imread(fullfile(d2(t).folder, d2(t).name)));
    a3 = im2double(imread(fullfile(d3(t).folder, d3(t).name)));
    a4 = im2double(imread(fullfile(d4(t).folder, d4(t).name)));

    I = (a1+a2+a3+a4)/4;

    % change ROI depending on oocyte posn.
    % I = imcrop(I,[5 5 1900 1900]);
    S = I;  % keep a copy for overlay display

    % --- Segmentation via contour_retardance pipeline (segment_oocyte) ---
    [BW_final, xc, yc, R_fit, poly, segCache, geom] = segment_oocyte(I, segParams, segCache);
    if isempty(poly); continue; end

    % Save mask
    BW_seq{t} = BW_final;
    % if mod(t,10)==0 && makeVideo
    %     imwrite(uint8(BW_final)*255, fullfile(outDir, sprintf('mask_%04d.png', t)));
    % end
    polySeq{t} = poly;
    geomSeq{t} = geom;
    xc_px(t) = xc;
    yc_px(t) = yc;
    R_fit_px(t) = R_fit;

    if isfield(geom, 'curvature_pxInv') && numel(geom.curvature_pxInv) == nBoundaryPts
        curvatureByArc_umInv(t,:) = geom.curvature_pxInv(:)' * px_per_um;
        radiusCurvatureByArc_um(t,:) = geom.radiusCurvature_px(:)' / px_per_um;
        thetaByArc_deg(t,:) = rad2deg(geom.theta_rad(:)');
        meanRadiusCurvature_um(t) = geom.meanRadiusCurvature_px / px_per_um;
        medianRadiusCurvature_um(t) = geom.medianRadiusCurvature_px / px_per_um;
    end

    % ---------- Map entire contour area to PIV nodes ----------
    % ---------- Map entire contour area to PIV nodes ----------
    X = Xc{t}; Y = Yc{t};

    % Erode mask so PIV nodes near the edge are excluded
    BW_eroded = imerode(BW_final, strel('disk', maskErodePx));
    B_eroded  = bwboundaries(BW_eroded);
    if ~isempty(B_eroded)
        [~, iLong] = max(cellfun(@(p) size(p,1), B_eroded));
        bE = B_eroded{iLong};
        polyMask = [bE(:,2), bE(:,1)];   % [x, y]
        if ~isequal(polyMask(1,:), polyMask(end,:)), polyMask(end+1,:) = polyMask(1,:); end
        INgrid = reshape(inpolygon(X(:), Y(:), polyMask(:,1), polyMask(:,2)), size(X));
    else
        INgrid = false(size(X));
    end
    INgridSeq{t} = INgrid;

    Ut = Uc{t}; Vt = Vc{t};
    Ut(~INgrid) = NaN;  Vt(~INgrid) = NaN;
    U_masked{t} = Ut;   V_masked{t} = Vt;

        % --- RAW IMAGE + PIV QUIVER OVERLAY (clean: no axes / title / text) ---
    [H,W] = size(S);

    % quiver node selection (guard against empty idx)
    idx = INgrid & ~isnan(Ut) & ~isnan(Vt);
    if ~any(idx(:))
        warning('Frame %d: no masked vectors to plot; showing all as fallback.', t);
        idx = ~isnan(Ut) & ~isnan(Vt);
    end

    % reduce arrow density: keep every arrowStride-th node in each direction
    keep = false(size(idx));
    keep(1:arrowStride:end, 1:arrowStride:end) = true;
    idx = idx & keep;

    % ---- ONE scale shared by the flow arrows AND the reference bar ----
    % a real flow arrow of scaleBar_refVel_umin will be exactly refLen long
    refVel_px = (scaleBar_refVel_umin/60) * dt * px_per_um;   % um/min -> px/frame
    refLen    = arrowScale * refVel_px;                        % bar length (px)

    if mod(t,10)==0
        figure('Name',sprintf('Overlay Frame %d',t), 'Color', 'w', 'Visible', 'off',...
            'Units','pixels','Position',[100 100 W H]);
        imagesc(S); 
        colormap gray;
        axis image; 
        axis off;
        % set(gca,'YDir','reverse'); 
        hold on;

        % flow vectors -- MANUAL scale (NOT autoscale) so they match the bar
        quiver(X(idx), Y(idx), arrowScale*Ut(idx), arrowScale*Vt(idx), 0, ...
            'Color', quiverColor, 'LineWidth', 2, 'MaxHeadSize', arrowHeadSize);

        % scale bar, anchored to the image's top-right corner
        xl = xlim;  yl = ylim;
        x2   = xl(2) - 0.15*diff(xl);
        x1   = x2 - refLen;
        yref = yl(1) + 0.07*diff(yl);

        if scaleBar_backdrop
            pX = 0.015*diff(xl);  pYt = 0.04*diff(yl);  pYb = 0.075*diff(yl);
            % patch([x1-pX, x2+pX, x2+pX, x1-pX], ...
            %       [yref-pYt, yref-pYt, yref+pYb, yref+pYb], ...
            %       'k', 'FaceAlpha',0.35, 'EdgeColor','none');
        end

        quiver(x1, yref, refLen, 0, 0, ...
            'Color', quiverColor, 'LineWidth', 2.5, 'MaxHeadSize', arrowHeadSize);
        text((x1+x2)/2, yref+0.03*diff(yl), sprintf('%.3g \\mum/min', scaleBar_refVel_umin), ...
            'Color', quiverColor, 'FontSize', 14, 'FontWeight','bold', ...
            'HorizontalAlignment','center', 'VerticalAlignment','top');

        exportgraphics(gca, fullfile(OvlerayDir, sprintf('PIV_and_AC_%04d.tif', t)));
        close(gcf);
    end


end

%% =================== Unit conversion to um/s =============================
um_per_px = 1/px_per_um;    % 0.16 um/px
for t = 1:nFrames
    U_masked{t} = U_masked{t} * um_per_px / dt;
    V_masked{t} = V_masked{t} * um_per_px / dt;
end
dx_um = dx_px * um_per_px;
dy_um = dy_px * um_per_px;

%% =================== Flow stats & divergence vs minutes ==================
vel_mean = nan(nFrames,1); vel_median = vel_mean; vel95 = vel_mean;
divRMS   = vel_mean;       divMed     = vel_mean;

for t = 1:nFrames
    U = U_masked{t}; V = V_masked{t}; IN = INgridSeq{t};
    Sp = hypot(U,V);
    m  = IN & ~isnan(Sp);
    if any(m(:))
        vals = Sp(m);
        vel_mean(t)   = mean(vals);
        vel_median(t) = median(vals);
        vel95(t)      = prctile(vals,95);
    end

    Uf = fillmissing(U,'nearest'); Vf = fillmissing(V,'nearest');
    [dUdx,~] = gradient(Uf, dx_um, dy_um);
    [~,dVdy] = gradient(Vf, dx_um, dy_um);
    dv = dUdx + dVdy;
    dvals = dv(IN);
    divRMS(t) = sqrt(mean(dvals.^2,'omitnan'));
    divMed(t) = median(dvals,'omitnan');
end



% Plots
figure('Color','k','Name','Mean Flow Speed <V(t)> ');
plot(t_min, vel_mean,'r','LineWidth',1.4); hold on;
xlabel('Time (min)'); ylabel('Mean Speed (um/s)'); legend('Mean','Location','best'); grid on;
title("Mean flow speed over time");
exportgraphics(gca, fullfile(outDir,'velocity_vs_time.png'), 'Resolution', 150);

figure('Color','k','Name','Divergence of V(t))');
plot(t_min, divMed,'b','LineWidth',1.4); ylabel('Mean div (s^{-1})');
xlabel('Time (min)'); grid on;
ylabel('Divergence (s^-1)');
title("Divergence over time");
exportgraphics(gca, fullfile(outDir,'divergence_vs_time.png'), 'Resolution', 150);

% Save CSV of metrics
R_fit_um = R_fit_px / px_per_um;
T = table(t_min(:), vel_mean(:), vel_median(:), vel95(:), divRMS(:), divMed(:), ...
    R_fit_um(:), meanRadiusCurvature_um(:), medianRadiusCurvature_um(:), ...
    'VariableNames', {'time_min','vel_mean_um_s','vel_median_um_s','vel95_um_s', ...
    'divRMS_s_inv','divMed_s_inv','R_fit_um','mean_radius_curvature_um','median_radius_curvature_um'});
writetable(T, fullfile(outDir, 'flow_metrics.csv'));


% Save MAT with masks, polygons, and masked fields
save(fullfile(outDir,'results_activecontour.mat'), ...
    'BW_seq','polySeq','geomSeq','INgridSeq','U_masked','V_masked', ...
    'vel_mean','vel_median','vel95','divRMS','divMed','t_min', ...
    'xc_px','yc_px','R_fit_px','R_fit_um', ...
    'curvatureByArc_umInv','radiusCurvatureByArc_um','thetaByArc_deg', ...
    'meanRadiusCurvature_um','medianRadiusCurvature_um', ...
    'px_per_um','dt','dx_um','dy_um','segParams','-v7.3');

fprintf('Done. Outputs saved under: %s\n', outDir);


% profile viewer


% make plot of peak speed and box size: too large and too small of a box
% should be zero


% vary box size and recheck all calibrations



% sanity checks: make sure wave speed of 45 um / min is reached. wave time
% is 8 mins and distance traveled is pi * R with R being approx ~ 200 um
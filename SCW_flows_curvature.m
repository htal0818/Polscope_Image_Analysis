close all;

theta=linspace(0,2*pi,101);  % for curvature bins

%% ========================== USER INPUTS ==================================

% test this on frames where SCW is actively occuring

% find a maximum speed the flows hit

% profile on

% --- Output root folder ---
outDir = '/Users/hridaytalreja/Downloads/Sep_oocytes_polscope_good_data/Jan15_2026_PIV/';
mkdir(outDir);



% % % %


base_dir = '/Users/hridaytalreja/Downloads/Sep_oocytes_polscope_good_data/Pos0/';


s1 = strcat(base_dir,'*_State1*');
s2 = strcat(base_dir,'*_State2*');
s3 = strcat(base_dir,'*_State3*');
s4 = strcat(base_dir,'*_State4*');

d1 = dir(s1);
d2 = dir(s2);
d3 = dir(s3);
d4 = dir(s4);


% % %


px_per_um = 6.25;       % 6.25 px per um (40x objective) or 3.12 for 20x
dt        = 30;

% --- Quiver appearance ---
quiverAutoScale      = true;
quiverAutoScaleFact  = 1.2;   % increase if arrows look small
quiverColor          = [0.2 0.85 0.2];  % green

% --- Save overlays ---
makeVideo = true;

%% ================== PIVlab arrays from workspace =========================

Uc = u_original; Vc = v_original; Xc = x; Yc = y;

nFrames = numel(Uc{1,:});

% grid spacing (pixels) for divergence after unit conversion
dx_px = median(diff(unique(Xc{1}(:))));
dy_px = median(diff(unique(Yc{1}(:))));


%% =================== Segmentation params (contour_retardance pipeline) ===
segParams.sigmaBlur            = 1;       % Gaussian blur sigma (px)
segParams.closeRadius          = 1;       % morphological close disk radius (px)
segParams.minArea              = 5000;    % min object area (px^2)
segParams.thresholdMode        = 'adaptive';  % 'adaptive', 'edge', or 'gradient'
segParams.segFromMask          = true;    % true for 4-state (oocyte is dark)
segParams.adaptSensitivity     = 0.7;
segParams.adaptNeighborhood    = 201;
segParams.edgeMethod           = 'Sobel';
segParams.edgeDilateRadius     = 2;
segParams.gradientPercentile   = 70;
segParams.useCaching           = true;
segParams.cacheIntensityThreshold = 0.02;
segParams.cacheForceRecalcEveryN  = 25;
segCache = [];

% --- Curvature computation parameters ---
polyOrder       = 50;            % polyfit order for r(theta) (matches kymograph)
nBoundary       = 500;           % samples along boundary for smooth polygon

%% =================== Containers =========================================
BW_seq     = cell(nFrames,1);
polySeq    = cell(nFrames,1);
INgridSeq  = cell(nFrames,1);
U_masked   = cell(nFrames,1);
V_masked   = cell(nFrames,1);
RADIUS_OF_CURVATURE = nan(nFrames, numel(theta)-1);
xc_saved = nan(nFrames,1);
yc_saved = nan(nFrames,1);

%% =================== MAIN LOOP: build mask (State-1) + map to PIV =======
for t = 1:nFrames
% for t = 75:77

    maskFile = fullfile(outDir, sprintf('mask_%04d.png', t));
    overlayFile = fullfile(outDir, sprintf('overlay_%04d.png', t));
    plotFile = fullfile(outDir, 'velocity_vs_time.png');
    csvFile  = fullfile(outDir, 'flow_metrics.csv');

    a1  = imread(fullfile(d1(t).folder, d1(t).name));
    a2 = imread(fullfile(d2(t).folder, d2(t).name));
    a3 = imread(fullfile(d3(t).folder, d3(t).name));
    a4 = imread(fullfile(d4(t).folder, d4(t).name));

    I = double(a1+a2+a3+a4);

    % change ROI depending on oocyte posn.
    I = imcrop(I,[5 5 1900 1900]);
    % Keep a copy for overlays
    Iraw = I;
    [H,W] = size(Iraw);

    % --- Segmentation via contour_retardance pipeline (segment_oocyte) ---
    [BW_final, xc, yc, R_fit, bndPoly, segCache] = segment_oocyte(I, segParams, segCache);
    if isempty(bndPoly); continue; end

    % --- Compute curvature from boundary (post-hoc polar polynomial fit) ---
    [RADIUS, poly] = compute_curvature_from_boundary(bndPoly, xc, yc, theta, polyOrder, nBoundary);

    % Clamp reconstructed boundary to image bounds and rebuild mask
    poly(:,1) = min(max(poly(:,1), 1), W);
    poly(:,2) = min(max(poly(:,2), 1), H);
    BW_final = poly2mask(poly(:,1), poly(:,2), H, W);
    BW_final = imfill(BW_final, 'holes');

    % Save mask and curvature
    BW_seq{t} = BW_final;
    RADIUS_OF_CURVATURE(t,:) = RADIUS;
    xc_saved(t) = xc;
    yc_saved(t) = yc;
    imwrite(uint8(BW_final)*255, fullfile(outDir, sprintf('mask_%04d.png', t)));

    polySeq{t} = poly;

    % ---------- Map entire contour area to PIV nodes ----------
    X = Xc{t}; Y = Yc{t};
    if ~isempty(poly)
        polyC = poly;
        if ~isequal(polyC(1,:), polyC(end,:)), polyC(end+1,:) = polyC(1,:); end
        INgrid = reshape(inpolygon(X(:), Y(:), polyC(:,1), polyC(:,2)), size(X));
    else
        INgrid = false(size(X));   % blank if we couldn't get a polygon
    end
    INgridSeq{t} = INgrid;

    Ut = Uc{t}; Vt = Vc{t};
    Ut(~INgrid) = NaN;  Vt(~INgrid) = NaN;
    U_masked{t} = Ut;   V_masked{t} = Vt;

    % --- QUIVER + MASK OVERLAY (aligned) ---
    S = Iraw;
    [H,W] = size(S);
    xMin = min(X(:)); xMax = max(X(:));
    yMin = min(Y(:)); yMax = max(Y(:));

    if mod(t,25)==0
        figure('Name',sprintf('Overlay Frame %d',t),'Color','w');
        imagesc([xMin xMax],[yMin yMax], S); axis image; colormap gray;
        set(gca,'YDir','reverse'); hold on;
    end

    % mask outline (converted to X,Y axes)
    Bmask = bwboundaries(BW_final);
    if ~isempty(Bmask)
        [~,ii] = max(cellfun(@(p) size(p,1), Bmask));
        b = Bmask{ii};
        sx = (xMax - xMin) / max(1,(W-1));
        sy = (yMax - yMin) / max(1,(H-1));
        bx = xMin + (b(:,2)-1)*sx;   % col -> x
        by = yMin + (b(:,1)-1)*sy;   % row -> y
        plot(bx, by, 'y-', 'LineWidth', 1.4);
    end

    % quiver (guard against empty idx)
    idx = INgrid & ~isnan(Ut) & ~isnan(Vt);
    if ~any(idx(:))
        warning('Frame %d: no masked vectors to plot; showing all as fallback.', t);
        idx = ~isnan(Ut) & ~isnan(Vt);
    end

    if exist('quiverAutoScaleFact','var') && ~isempty(quiverAutoScaleFact)
        quiver(X(idx), Y(idx), Ut(idx), Vt(idx), ...
            'AutoScale','on','AutoScaleFactor', quiverAutoScaleFact, ...
            'Color', quiverColor, 'LineWidth', 1.2, 'MaxHeadSize', 1.1);
    else
        quiver(X(idx), Y(idx), Ut(idx), Vt(idx), 0, ...
            'Color', quiverColor, 'LineWidth', 1.2, 'MaxHeadSize', 1.1);
    end

    text(xMin+10, yMin+20, sprintf('Frame %d', t), ...
        'Color','w','FontSize',14,'FontWeight','bold');

    title(sprintf('Quiver + Boundary Estimation — Frame %d', t));

    if mod(t,25)==0
        exportgraphics(gca, fullfile(outDir, sprintf('PIV_map_%04d.png', t)), 'Resolution', 150);
        drawnow;
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

t_sec = (0:nFrames-1) * dt;
t_min = t_sec/60;

% Plots
figure('Color','w','Name','Mean Flow Speed <V(t)> ');
plot(t_min, vel_mean,'--','LineWidth',1.4); hold on;
xlabel('Time (min)'); ylabel('Mean Speed (um/s)'); legend('Mean','Location','best'); grid on;
exportgraphics(gca, fullfile(outDir,'velocity_vs_time.png'), 'Resolution', 150);

figure('Color','w','Name','Divergence of V(t))');
plot(t_min, divMed,'--','LineWidth',1.4); ylabel('Mean div (s^{-1})');
xlabel('Time (min)'); grid on;
ylabel('Divergence (s^-1)');
exportgraphics(gca, fullfile(outDir,'divergence_vs_time.png'), 'Resolution', 150);

% Save CSV of metrics
T = table(t_min(:), vel_mean(:), vel_median(:), vel95(:), divRMS(:), divMed(:), ...
    'VariableNames', {'time_min','vel_mean_um_s','vel_median_um_s','vel95_um_s','divRMS_s_inv','divMed_s_inv'});
writetable(T, fullfile(outDir, 'flow_metrics.csv'));


% Save MAT with masks, polygons, and masked fields
save(fullfile(outDir,'results_curvature.mat'), ...
'BW_seq','polySeq','INgridSeq','U_masked','V_masked', ...
'vel_mean','vel_median','vel95','divRMS','divMed','t_min', ...
'RADIUS_OF_CURVATURE','xc_saved','yc_saved', ...
'px_per_um','dt','dx_um','dy_um','-v7.3');

fprintf('Done. Outputs saved under: %s\n', outDir);


% profile viewer


% make plot of peak speed and box size: too large and too small of a box
% should be zero


% vary box size and recheck all calibrations



% sanity checks: make sure wave speed of 45 um / min is reached. wave time
% is 8 mins and distance traveled is pi * R with R being approx ~ 200 um

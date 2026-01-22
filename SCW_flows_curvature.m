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


px_per_um = 6.25;       % 6.25 px per µm (40x objective) or 3.12 for 20x
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





%% =================== Mask params =========================
maskMethod      = 'curvature';   % 'curvature' (recommended) or 'activecontour'
sigmaBlur       = 1.0;           % pre-blur (pixels); keep small to preserve edge
threshFrac      = 0.85;          % mask = I < threshFrac*mean2(I)
se              = strel('diamond',5);

polyOrder       = 50;            % polyfit order for r(theta) (matches kymograph)
nBoundary       = 500;           % samples along boundary for smooth polygon
pixels_to_average = 0;           % optional radial averaging for sampling (not needed for mask)

% If you ever want to fall back to activecontour, keep these here (unused unless maskMethod='activecontour')
acItersMax      = 300;
acChunk         = 50;
acSmooth        = 2;
acContraction   = 0.25;
annulusMin      = 0.15;
annulusMax      = 0.65;

%% =================== Containers =========================================
BW_seq     = cell(nFrames,1);
polySeq    = cell(nFrames,1);
INgridSeq  = cell(nFrames,1);
U_masked   = cell(nFrames,1);
V_masked   = cell(nFrames,1);
RADIUS_OF_CURVATURE = nan(nFrames, numel(theta)-1);
xc_saved = nan(nFrames,1);
yc_saved = nan(nFrames,1);


prevBW = []; prevC = [];

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

    I= ((a1+a2+a3+a4));

    xmin = 5;
    ymin = 5;
    width = 10;
    height = 10;

    

    
    

    % change ROI depdending on oocyte posn.
    I=imcrop(I,[5 5 1900 1900]);
    % Keep a copy for overlays
    Iraw = I;

    % --- Curvature-based boundary reconstruction (copied as closely as possible from kymograph.m) ---
    [H,W] = size(Iraw);

    if strcmpi(maskMethod,'curvature')

        % 1) coarse threshold mask (same logic as kymograph.m)
        % NOTE: kymograph.m computes a gaussian-blurred image but thresholds the *raw* summed state image.


        Iraw = imgaussfilt(Iraw,sigmaBlur);


        BW0 = Iraw < (threshFrac*mean2(Iraw));

        BW0 = imdilate(BW0, se);
        BW0 = imfill(BW0,'holes');
        BW0 = imerode(BW0, se);
       

        % keep the largest connected component (kymograph uses bwlabel + regionprops Area loop)
        labels = bwlabel(BW0);
        clear Area
        for i =1:max(max(labels))
            temp=regionprops(labels==i,'Area');
            Area(i)=temp.Area;
        end

        mask=labels==find(Area==max(max(Area)));

        % 2) edge pixels of the coarse mask (kymograph: edges=imgradient(mask); edges=edges>0)
        edges = imgradient(BW0);
        edges = edges > 0;

        
        xx=[];yy=[];
        for i = 1:size(I,1)
            for j = 1:size(I,2)
                if edges(i,j)==1
                    xx=[xx j];
                    yy=[yy i];
                end
            end
        end
        clear angle r

        % defaults
        RADIUS = nan(1, numel(theta)-1);
        poly   = [];
        xc = W/2; yc = H/2;

        if numel(xx) < 50
            warning('Frame %d: too few edge pixels for circfit; using BW0.', t);
            BW_final = BW0;
        else
            % --- circle fit to get center (kymograph: [xc yc R] = circfit(x,y))
            [xc yc R] = circfit(xx,yy);

            xc_saved(t) = xc;
            yc_saved(t) = yc;

            % --- compute r and angle for each edge point (kymograph logic with acos + y>yc branch)
            nPts = numel(xx);
            r = zeros(1,nPts);
            angle = zeros(1,nPts);

            for kk = 1:nPts
                r(kk) = norm([xx(kk)-xc, yy(kk)-yc]);

                if yy(kk) > yc
                    angle(kk) = acos( dot( ([xx(kk)-xc, yy(kk)-yc]) / norm([xx(kk)-xc, yy(kk)-yc]), [1 0]) );
                else
                    angle(kk) = 2*pi - acos( dot( ([xx(kk)-xc, yy(kk)-yc]) / norm([xx(kk)-xc, yy(kk)-yc]), [1 0]) );
                end
            end

            % --- FIRST PASS: fit r(angle) and fill the middle half of RR (kymograph indices 126:375)
            RRrow = nan(1, nBoundary);

            param = polyfit(angle, r, polyOrder);   % kymograph: polyfit(angle,r,50)
            xgrid = linspace(0, 2*pi, nBoundary);   % kymograph: x=linspace(0,2*pi,500)
            y1    = polyval(param, xgrid);

            r_theta_p   = polyder(param);
            r_2theta_p  = polyder(r_theta_p);
            r_theta     = polyval(r_theta_p, xgrid);
            r_2theta    = polyval(r_2theta_p, xgrid);

            % kymograph fills curvature bins only for i = 25:75 (for 101 theta)
            for ii = (length(theta)-1)/2 - (length(theta)-1)/4 : (length(theta)-1)/2 + (length(theta)-1)/4
                sel = xgrid > theta(ii) & xgrid < theta(ii+1);
                if any(sel)
                    num = ((y1(sel).^2 + r_theta(sel).^2).^(3/2));
                    den = abs(y1(sel).^2 + 2*r_theta(sel).^2 - y1(sel).*r_2theta(sel));
                    RADIUS(ii) = mean(num ./ max(den, eps));
                end
            end

            RRrow(126:375) = y1(126:375);

            % --- SECOND PASS: sort and shift by pi to handle wrap-around (kymograph logic)
            [angle2, Iord] = sort(angle);
            r2 = r(Iord);

            angle2 = angle2 + pi;
            angle2(angle2 > 2*pi) = angle2(angle2 > 2*pi) - 2*pi;

            param2 = polyfit(angle2, r2, polyOrder);
            x2     = linspace(0, 2*pi, nBoundary);
            y2     = polyval(param2, x2);

            r_theta_p2  = polyder(param2);
            r_2theta_p2 = polyder(r_theta_p2);
            r_theta2    = polyval(r_theta_p2, x2);
            r_2theta2   = polyval(r_2theta_p2, x2);

            x2 = x2 - pi;
            x2(x2 < 0) = 2*pi + x2(x2 < 0);

            y2 = circshift(y2, nBoundary/2);   % kymograph: circshift(y1,250) for 500 samples

            RRrow(1:125)   = y2(1:125);
            RRrow(376:500) = y2(376:500);

            % curvature bins for the remaining quarters (kymograph does this in two loops)
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

            % --- Reconstruct final smooth boundary from RR (kymograph: XX=RR.*cos(x)+xc, YY=RR.*sin(x)+yc)
            xfinal = linspace(0, 2*pi, nBoundary);
            XX = RRrow .* cos(xfinal) + xc;
            YY = RRrow .* sin(xfinal) + yc;

            % clamp to image bounds for poly2mask stability
            XX = min(max(XX,1), W);
            YY = min(max(YY,1), H);

            % 3) build the final mask from the reconstructed smooth boundary
            BW_final = poly2mask(XX, YY, H, W);
            BW_final = imfill(BW_final,'holes');

            poly = [XX(:) YY(:)];
        end

    else
        error('maskMethod="%s" not implemented in this version. Set maskMethod="curvature".', maskMethod);
    end

% Save mask
    BW_seq{t} = BW_final;
    RADIUS_OF_CURVATURE(t,:) = RADIUS;
    xc_saved(t) = xc;
    yc_saved(t) = yc;
    imwrite(uint8(BW_final)*255, fullfile(outDir, sprintf('mask_%04d.png', t)));
    % imshow(BW_final)
    
  


    


    % Polygon for mapping/outline (use curvature polygon if available)
    if isempty(poly)
        B = bwboundaries(BW_final);
        if ~isempty(B)
            [~,ii] = max(cellfun(@(p) size(p,1), B));
            b   = B{ii}; step = max(1, floor(size(b,1)/400));
            bds = b(1:step:end, :);
            poly = [bds(:,2), bds(:,1)];
        else
            poly = [];
        end
    end
    polySeq{t} = poly;


    % ---------- Map entire contour area to PIV nodes ----------
    X = Xc{t}; Y = Yc{t};
    if ~isempty(poly)
        if ~isequal(poly(1,:), poly(end,:)), poly(end+1,:) = poly(1,:); end
        INgrid = reshape(inpolygon(X(:), Y(:), poly(:,1), poly(:,2)), size(X));
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

%% =================== Unit conversion to µm/s =============================
um_per_px = 1/px_per_um;    % 0.16 µm/px
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

t_sec = (0:nFrames-1) * dt;  % 15 s = 0.25 min
t_min = t_sec/60;



% Plots
figure('Color','w','Name','Mean Flow Speed <V(t)> ');
plot(t_min, vel_mean,'--','LineWidth',1.4); hold on;
xlabel('Time (min)'); ylabel('Mean Speed (µm/s)'); legend('Mean','Location','best'); grid on;
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
'px_per_um','dt','dx_um','dy_um','-v7.3');

fprintf('Done. Outputs saved under: %s\n', outDir);


% profile viewer


% make plot of peak speed and box size: too large and too small of a box
% should be zero


% vary box size and recheck all calibrations 



% sanity checks: make sure wave speed of 45 um / min is reached. wave time
% is 8 mins and distance traveled is pi * R with R being approx ~ 200 um


%% =================== Local helper: circle fit (no toolbox dependency) ====
% function [xc,yc,R] = circfit(x,y)
% % CIRCfit: simple algebraic circle fit (Kasa method).
% % Inputs x,y are column or row vectors of edge pixel coordinates.
% x = x(:); y = y(:);
% A = [2*x, 2*y, ones(size(x))];
% b = x.^2 + y.^2;
% sol = A\b;
% xc = sol(1);
% yc = sol(2);
% c  = sol(3);
% R  = sqrt(max(c + xc^2 + yc^2, 0));
% end


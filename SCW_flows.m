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

state1_input = '/Users/hridaytalreja/Desktop/Jan_data_2026/2026_01_07_test_20x_FSW_and_alignment/SMS_2026_0107_1425_1/state1/';   

% --- Output root folder ---
outDir = '/Users/hridaytalreja/Desktop/Jan_data_2026/2026_01_07_test_20x_FSW_and_alignment/SMS_2026_0107_1425_1/PIV_Jan13_2026/';
mkdir(outDir);
% mkdir(outDir);


% % % %


base_dir = '/Users/hridaytalreja/Desktop/Jan_data_2026/2026_01_07_test_20x_FSW_and_alignment/SMS_2026_0107_1425_1/Pos0/';


s1 = strcat(base_dir,'/*State1*');
s2 = strcat(base_dir,'/*State2*');
s3 = strcat(base_dir,'/*State3*');
s4 = strcat(base_dir,'/*State4*');

d1 = dir(s1);
d2 = dir(s2);
d3 = dir(s3);
d4 = dir(s4);


% % % 


px_per_um = 6.25/2;       % 6.25 px per µm (40x objective) or 3.12 for 20x
dt        = 15;         % 30 s per frame

% --- Quiver appearance ---
quiverAutoScale      = true;
quiverAutoScaleFact  = 1.2;   % increase if arrows look small
quiverColor          = [0.2 0.85 0.2];  % green

% --- Save overlays ---
makeVideo = true;     

%% ================== PIVlab arrays from workspace =========================

Uc = u_original; Vc = v_original; Xc = x; Yc = y;

nFrames = numel(Uc);

% grid spacing (pixels) for divergence after unit conversion
dx_px = median(diff(unique(Xc{1}(:))));
dy_px = median(diff(unique(Yc{1}(:))));





%% =================== Active contour params / reuse / early stop ==========
sigmaBlur       = 3.0;     % pre-blur for stability
minAreaPx       = 200;     % min component size (pixels)
fudgeFactor     = 0.5;     % Sobel threshold scale
se_len          = 3;       % dilation line length
smooth_iters    = 2;       % diamond erosions at the end

acItersMax      = 300;     % max iterations cap
acChunk         = 50;      % iterate in chunks (check convergence)
acSmooth        = 2;       % smoothing 
acContraction   = 0.25;       % >0 inward bias, <0 outward, 0 neutral

growPx          = 0;       % dilate previous mask for growth
shiftWithCentroid = true;  % translate prev mask toward current seed

iouStop         = 0.995;   % early stop if IoU ≥ this
fracChangeStop  = 1e-4;    % or fractional pixel change ≤ this

annulusMin      = 0.15;    % ROI inner radius (fraction of min(H,W))
annulusMax      = 0.65;    % ROI outer radius (fraction of min(H,W))

% se90 = strel('line', se_len, 90);
% se0  = strel('line', se_len, 0);
se  = strel('diamond', 5);

%% =================== Containers =========================================
BW_seq     = cell(nFrames,1);
polySeq    = cell(nFrames,1);
INgridSeq  = cell(nFrames,1);
U_masked   = cell(nFrames,1);
V_masked   = cell(nFrames,1);

prevBW = []; prevC = [];

%% =================== MAIN LOOP: build mask (State-1) + map to PIV =======
for t = 150:nFrames
% for t = 75:77



    maskFile = fullfile(outDir, sprintf('mask_%04d.png', t));
    overlayFile = fullfile(outDir, sprintf('overlay_%04d.png', t));
    plotFile = fullfile(outDir, 'velocity_vs_time.png');
    csvFile  = fullfile(outDir, 'flow_metrics.csv');

 


    a1  = im2double(imread(fullfile(d1(t).folder, d1(t).name)));


    a2 = im2double(imread(fullfile(d2(t).folder, d2(t).name)));
    a3 = im2double(imread(fullfile(d3(t).folder, d3(t).name)));
    a4 = im2double(imread(fullfile(d4(t).folder, d4(t).name)));

    xmin = 5;
    ymin = 5;
    width = 10;
    height = 10;


    I= (a1+a2+a3+a4)/1;
    % a = im2gray(a);

    % change ROI depdending on oocyte posn.
    I=imcrop(I,[5 5 1900 1900]);
  

    % --- 1) Blur ---
    I = imgaussfilt(I, sigmaBlur);

    % --- 2) Gradient magnitude & threshold ---

    mask = I<.85*mean2(I);
    mask = imdilate(mask,se);
    mask=imfill(mask,'holes');
    edges = imgradient(mask);
    % gThr      = prctile(Gmag(:), 85);
    edges=edges>0;

    % Clean up init (keep strong rim)
    % BW_init = imfill(BW_grad, 'holes');
    % BW_init = bwareaopen(BW_init, 50);
    % if any(BW_init(:)), BW_init = bwareafilt(BW_init, 1); end

    % Smooth/close, then dilate

    % BW_init = imclose(BW_init, strel('disk', 3, 0));

    BW_init = edges;
    % --- centroid from current seed (or center fallback) ---
    if any(BW_init(:))
        Sreg = regionprops(BW_init, 'Area','Centroid');
        [~,imax] = max([Sreg.Area]);
        Cseed = Sreg(imax).Centroid;   % [x y]
    else
        Cseed = [size(I,2) size(I,1)]/2;  % [x y]
    end

    % --- ROI annulus to keep AC near the oocyte ---
    [H,W] = size(I);
    if ~isempty(prevC), C = prevC; else, C = Cseed; end
    cx = C(1); cy = C(2);
    [Xg,Yg] = meshgrid(1:W,1:H);
    r = sqrt((Xg-cx).^2 + (Yg-cy).^2);
    rmin = annulusMin*min(H,W);  rmax = annulusMax*min(H,W);
    ROI  = (r >= rmin) & (r <= rmax);

    % --- reuse previous frame's mask as seed ---
    BW_seed = BW_init;
    if ~isempty(prevBW)
        BW_prev = prevBW;
        if shiftWithCentroid && ~isempty(prevC)
            dxy = round(Cseed - prevC);  % [dx dy]
            if any(dxy~=0)
                try
                    BW_prev = imtranslate(BW_prev, dxy, 'FillValues', 0, 'OutputView', 'same');
                catch
                    % older MATLAB: ignore shift if not supported
                end
            end
        end
        BW_prev = imdilate(BW_prev, strel('disk', growPx));  % allow modest growth
        BW_seed = (BW_seed | (BW_prev & ROI));                % fuse & bind to ROI
    end
    if ~any(BW_seed(:))
        BW_seed(round(H/2), round(W/2)) = true;
        BW_seed = imdilate(BW_seed, strel('disk', 5));
    end



    BW_ac0  = BW_init;

    BW_ac0 = imdilate(BW_init, strel('disk', 2, 0));


    % --- 3) Active contour refinement ---

    % --- finalize mask; keep largest component; remember for next frame ---
    BW_ac = bwareafilt(BW_ac0, 1);
    BW_ac = imfill(BW_ac, 'holes');
    BW_ac = imopen(BW_ac, strel('disk', 2, 0));

    if any(BW_ac(:))
        Sreg = regionprops(BW_ac, 'Area','Centroid');
        [~,imax] = max([Sreg.Area]);
        BW_final = BW_ac;   % already single blob
        C_final  = Sreg(imax).Centroid;
        prevBW   = BW_final;   % <— reuse next frame
        prevC    = C_final;    % <— reuse next frame
    else
        if ~isempty(prevBW)
            BW_final = prevBW;   % fallback to last good
        else
            BW_final = false(size(I));
        end
    end

    I_for_ac = imgaussfilt(I, sigmaBlur);
    I_for_ac = I_for_ac.*ROI + (~ROI).*median(I_for_ac(:));

    BW_ac = BW_seed;
    BW_prev_iter = BW_ac;
    nChunks = ceil(acItersMax / acChunk);

    for c = 1:nChunks
        try
            BW_ac = activecontour(I_for_ac, BW_ac, acChunk, 'edge', ...
                'SmoothFactor', acSmooth, 'ContractionBias', acContraction);
        catch
            BW_ac = activecontour(I_for_ac, BW_ac, acChunk, 'edge');
        end
        BW_ac = BW_ac & ROI;
        BW_ac = imfill(BW_ac,'holes');
        BW_ac = imopen(BW_ac, strel('disk',1));

        inter = nnz(BW_ac & BW_prev_iter);
        uni   = nnz(BW_ac | BW_prev_iter);
        iou   = inter / max(1, uni);
        fracChange = nnz(xor(BW_ac, BW_prev_iter)) / numel(BW_ac);

        if (iou >= iouStop) || (fracChange <= fracChangeStop)
            break  % early convergence
        end
        BW_prev_iter = BW_ac;
    end



    % Final clean-up
    BW_final = bwareafilt(BW_ac, 1);
    BW_final = imfill(BW_final, 'holes');
    BW_final = imopen(BW_final, strel('disk', 2, 0));

    % Save mask
    BW_seq{t} = BW_final;
    imwrite(uint8(BW_final)*255, fullfile(outDir, sprintf('mask_%04d.png', t)));
    % imshow(BW_final)
    
  


    


    % Polygon for mapping/outline
    B = bwboundaries(BW_final);
    if ~isempty(B)
        [~,ii] = max(cellfun(@(p) size(p,1), B));
        b   = B{ii}; step = max(1, floor(size(b,1)/400));
        bds = b(1:step:end, :);
        poly = [bds(:,2), bds(:,1)];
    else
        poly = [];
    end
    polySeq{t} = poly;

    % ---------- Map entire contour area to PIV nodes ----------
    X = Xc{t}; Y = Yc{t};
    if ~isempty(poly)
        if ~isequal(poly(1,:), poly(end,:)), poly(end+1,:) = poly(1,:); end
        INgrid = reshape(inpolygon(X(:), Y(:), poly(:,1), poly(:,2)), size(X));
    else
        INgrid = false(size(X));   % blank if we couldn’t get a polygon
    end
    INgridSeq{t} = INgrid;

    Ut = Uc{t}; Vt = Vc{t};
    Ut(~INgrid) = NaN;  Vt(~INgrid) = NaN;
    U_masked{t} = Ut;   V_masked{t} = Vt;

    % --- QUIVER + MASK OVERLAY (aligned) ---
    [H,W] = size(S);
    xMin = min(X(:)); xMax = max(X(:));
    yMin = min(Y(:)); yMax = max(Y(:));

    
    
    if mod(t,10)==0
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

    title(sprintf('Quiver + Active Contour Mask — Frame %d', t));

    if mod(t,10)==0
        exportgraphics(gca, fullfile(outDir, sprintf('PIV_and_AC_%04d.png', t)), 'Resolution', 150);
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

t_min = (0:nFrames-1)*0.25;  % 15 s = 0.25 min

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
save(fullfile(outDir,'results_activecontour.mat'), ...
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



close all; clear all; clc


function [K, angles_deg, time_s] = minimal_polscope_bins_final(inputPath, outFolder, opts)
% Minimal PolScope pipeline (plots only after the loop).
% - inputPath: multipage .tif OR a folder of .tif images
% - outFolder: where to save outputs
% - opts (optional): struct with fields:
%     .Ntheta (default 100)         number of angular bins over 0..360°
%     .frameInterval (default 30)   seconds per frame
%     .sigmaBlur (default 2)        Gaussian sigma (px)
%     .minArea (default 200)        min object area (px^2)
%     .pxPerMicron (default 6.25)   pixels per micron (=> 0.16 µm/px)
%     .saveOverlays (default false) save overlay PNGs per frame (headless)
%
% Outputs:
%   K           [T × Ntheta] kymograph (boundary intensity vs angle)
%   angles_deg  [1 × Ntheta] angle axis (degrees)
%   time_s      [T × 1]      time axis (seconds)

% ---------- defaults ----------
if nargin < 3, opts = struct; end
function val = getOpt(opts, field, default)
    if isfield(opts, field) && ~isempty(opts.(field))
        val = opts.(field);
    else
        val = default;
    end
end

Ntheta        = getOpt(opts,'Ntheta', 100);
frameInterval = getOpt(opts,'frameInterval', 30);
sigmaBlur     = getOpt(opts,'sigmaBlur', 2);
minArea       = getOpt(opts,'minArea', 200);
pxPerMicron   = getOpt(opts,'pxPerMicron', 6.25);
saveOverlays  = getOpt(opts,'saveOverlays', false);


umPerPx = 1 / pxPerMicron;

if ~exist(outFolder,'dir'), mkdir(outFolder); end

% ---------- read files ----------
if isfolder(inputPath)
    d = dir(fullfile(inputPath,'*.tif'));
    if isempty(d), error('No .tif files found in folder: %s', inputPath); end
    [~,ord] = sort({d.name});
    d = d(ord);
    T = numel(d);
    I0 = im2double(imread(fullfile(d(1).folder, d(1).name)));
    [H,W] = size(I0);
    readFrame = @(t) im2double(imread(fullfile(d(t).folder, d(t).name)));
    frameName = @(t) d(t).name;
else
    info = imfinfo(inputPath);
    T = numel(info); if T < 1, error('No frames in: %s', inputPath); end
    H = info(1).Height; W = info(1).Width;
    readFrame = @(t) im2double(imread(inputPath, t));
    frameName = @(t) sprintf('%s (page %d)', inputPath, t);
end

% ---------- storage ----------
time_s     = (0:T-1)' * frameInterval;
angles_deg = linspace(0, 360, Ntheta);
K          = nan(T, Ntheta);


prevBW  = [];
prevC   = [W/2, H/2];
prevRow = nan(1, Ntheta);

% ---------- processing ----------
for t = 1:T
    I = readFrame(t);

    %  blur
    I_blur = imgaussfilt(I, sigmaBlur);

    % Otsu segmentation 
    Totsu = graythresh(I_blur);
    BW = I_blur > Totsu;
    BW = imfill(BW,'holes');
    BW = bwareaopen(BW, minArea);

    if ~any(BW(:))
        
        [Gmag,~] = imgradient(I_blur);
        thrG = max(2*mean(Gmag(:)), prctile(Gmag(:),80));
        BW   = Gmag >= thrG;
        BW   = imfill(BW,'holes');
        BW   = bwareaopen(BW, minArea);
    end

    % keep largest component and reuse previous BW if empty
    L = bwlabel(BW,8);
    if max(L(:)) >= 1
        S = regionprops(L, 'Area','Centroid');
        [~,iMax] = max([S.Area]);
        BW = (L == iMax);
        C  = S(iMax).Centroid;   % [x y]
    else
        BW = prevBW;
        C  = prevC;
    end

    
    if saveOverlays && any(BW(:))
        B = bwboundaries(BW);
        if ~isempty(B)
            [~,ii] = max(cellfun(@(p) size(p,1), B));
            boundary = B{ii};
            fig = figure('Visible','off'); imshow(I,[]); hold on; axis image
            plot(boundary(:,2), boundary(:,1), 'r-', 'LineWidth', 1.4);
            plot(C(1), C(2), 'g+', 'LineWidth', 1.2, 'MarkerSize', 10);
            title(sprintf('Overlay — %s', frameName(t)),'Interpreter','none');
            exportgraphics(gca, fullfile(outFolder, sprintf('overlay_%03d.png', t)));
            close(fig);
        end
    end

     
    row = nan(1, Ntheta);
    if any(BW(:))
        B = bwboundaries(BW);
        if ~isempty(B)
            [~,ii] = max(cellfun(@(p) size(p,1), B));     
            bnd = B{ii};                                  % [y, x]
            yb  = bnd(:,1); xb = bnd(:,2);

            th = atan2(yb - C(2), xb - C(1));
            th(th < 0) = th(th < 0) + 2*pi;               

            F  = griddedInterpolant({1:H,1:W}, I, 'linear','nearest');
            ib = F(yb, xb);                                % intensity on boundary

            edges = linspace(0, 2*pi, Ntheta+1);          
            bin   = discretize(th, edges);                
            valid = ~isnan(bin);

            if any(valid)
                sumBins   = accumarray(bin(valid), ib(valid), [Ntheta 1], @nansum, NaN);
                countBins = accumarray(bin(valid), 1,          [Ntheta 1], @nansum, NaN);
                row = (sumBins ./ countBins).';
            end
        end
    end


    if all(isnan(row))
        if any(BW(:))
            mval = mean(I(BW),'omitnan');
            if isfinite(mval), row(:) = mval; else, row = prevRow; end
        else
            row = prevRow;
        end
    else
        % circular interpolation for missing bins
        bad = isnan(row);
        if any(bad)
            goodIdx = find(~bad);
            if numel(goodIdx) >= 2
                xw = [goodIdx, goodIdx + Ntheta];
                vw = [row(goodIdx), row(goodIdx)];
                xi = 1:Ntheta;
                yi = interp1(xw, vw, xi, 'linear', 'extrap');
                m  = isnan(yi);
                if any(m), yi(m) = interp1(xw, vw, xi(m) + Ntheta, 'linear', 'extrap'); end
                row = yi;
            else
                mval = mean(I(BW),'omitnan');
                if isfinite(mval), row(:) = mval; end
            end
        end
    end





K(t,:) = row;

    % prevBW  = BW;
    % prevC   = C;
    % prevRow = row;
end

% ---------- PLOTS  ----------
% Kymograph: intensity vs time & angle
figK = figure('Visible','on');
imagesc(1:size(K,2), time_s, K); set(gca,'YDir','normal');
xlabel('Angle bin'); ylabel('Time (s)');
title('Retardance at cortex (time vs. angle)');
colormap parula; colorbar;
exportgraphics(gca, fullfile(outFolder,'kymograph_time_vs_angle_100bins.png'));
close(figK);

% Radial intensity vs distance 
frames_pick = unique([1, max(1,round(T/2)), T]);
dr_um = 0.5;  % radial sampling 

figR = figure('Visible','on'); hold on
for idx = 1:numel(frames_pick)
    t = frames_pick(idx);
    I = readFrame(t);

    
    I_blur = imgaussfilt(I, sigmaBlur);
    Totsu  = graythresh(I_blur);
    BWt    = imfill(bwareaopen(I_blur > Totsu, minArea), 'holes');
    if ~any(BWt(:))
        [Gmag,~] = imgradient(I_blur);
        thrG = max(2*mean(Gmag(:)), prctile(Gmag(:),85));
        BWt  = imfill(bwareaopen(Gmag >= thrG, minArea), 'holes');
    end
    L = bwlabel(BWt,8);
    if max(L(:)) >= 1
        S = regionprops(L,'Area','Centroid');
        [~,iMax] = max([S.Area]);
        BWt = (L == iMax);
        C   = S(iMax).Centroid;
    else
        C   = [W/2, H/2];
    end

    % determine a radius fully inside image
    Rmax_safe = floor(min([C(1)-1, W-C(1), C(2)-1, H-C(2)]));
    if Rmax_safe < 5, Rmax_safe = min([H W])/4; end
    R_um = 0:dr_um:Rmax_safe*umPerPx;
    R_px = R_um / umPerPx;

    % angle-average intensity at each radius
    theta = linspace(0, 2*pi, 360+1); theta(end) = [];
    [RR,TT] = ndgrid(R_px, theta);
    XX = C(1) + RR.*cos(TT);  YY = C(2) + RR.*sin(TT);
    Ipol = interp2(1:W, 1:H, I, XX, YY, 'linear', NaN);
    radial_mean = mean(Ipol, 2, 'omitnan');

    plot(R_um, radial_mean, 'LineWidth', 1.6, ...
        'DisplayName', sprintf('t = %.0f s', time_s(t)));
end
xlabel('Radius from center (\mum)');
ylabel('Retardance (nm)');
title('Angle-averaged radial intensity');
grid on; legend show;
exportgraphics(gca, fullfile(outFolder,'radial_intensity_selected_times.png'));
close(figR);

% ---------- save data ----------
save(fullfile(outFolder,'minimal_polscope_bins_final.mat'), ...
     'K','angles_deg','time_s','Ntheta','frameInterval','sigmaBlur','minArea','pxPerMicron');

fprintf('Done. Saved outputs to: %s\n', outFolder);
end



% Input call with a folder of TIF files
[K, ang_deg, tsec] = minimal_polscope_bins_final( ...
    '/Users/hridaytalreja/Downloads/retardance_ims/', ...
    '/Users/hridaytalreja/Downloads/retardance_ims/output_test', ...
    struct('Ntheta',100,'frameInterval',30,'sigmaBlur',2,'minArea',200,'saveOverlays',false));




%
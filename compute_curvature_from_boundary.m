function [RADIUS, poly, param, param2, r_theta_p, r_theta_p2, theta_boundary, RRrow] = ...
    compute_curvature_from_boundary(bndPoly, xc, yc, theta, polyOrder, nBoundary)
% COMPUTE_CURVATURE_FROM_BOUNDARY  Polar polynomial fit + curvature from boundary points.
%
% Takes boundary polygon points (e.g. from segment_oocyte / bwboundaries),
% converts to polar coordinates, fits two-pass polynomials to handle wrap-
% around, computes radius of curvature per angular bin, and reconstructs a
% smooth boundary polygon.
%
% This logic is extracted from the inline curvature code formerly in
% SCW_flows_curvature.m and boundary_flows.m (make_oocyte_mask_curvature).
%
% Inputs:
%   bndPoly    - [N x 2] boundary points [x, y] (e.g. from bwboundaries)
%   xc, yc     - circle fit center coordinates
%   theta      - angular bin edges, e.g. linspace(0, 2*pi, 101)
%   polyOrder  - polynomial order for r(theta) fit (default: 50)
%   nBoundary  - number of samples for smooth boundary (default: 500)
%
% Outputs:
%   RADIUS          - [1 x (numel(theta)-1)] radius of curvature per angular bin
%   poly            - [nBoundary x 2] reconstructed smooth boundary [XX, YY]
%   param           - first-pass polynomial coefficients (for tangent computation)
%   param2          - second-pass polynomial coefficients (for tangent computation)
%   r_theta_p       - first-pass derivative polynomial coefficients
%   r_theta_p2      - second-pass derivative polynomial coefficients
%   theta_boundary  - [1 x nBoundary] theta values at boundary points
%   RRrow           - [1 x nBoundary] fitted radius values
%
% REQUIREMENTS:
%   - circfit.m (for the calling scripts)

if nargin < 5 || isempty(polyOrder); polyOrder = 50; end
if nargin < 6 || isempty(nBoundary); nBoundary = 500; end

%% --- Default outputs ---
nThetaBins = numel(theta) - 1;
RADIUS = nan(1, nThetaBins);
poly = [];
param = []; param2 = [];
r_theta_p = []; r_theta_p2 = [];
theta_boundary = [];
RRrow = nan(1, nBoundary);

%% --- Extract boundary point coordinates ---
xx = bndPoly(:,1)';  % x coordinates (column -> row)
yy = bndPoly(:,2)';  % y coordinates (column -> row)

if numel(xx) < 50
    warning('compute_curvature_from_boundary: Too few boundary points (%d) for polynomial fit.', numel(xx));
    return;
end

%% --- Convert to polar coordinates relative to (xc, yc) ---
nPts = numel(xx);
r = zeros(1, nPts);
angle = zeros(1, nPts);

for kk = 1:nPts
    r(kk) = norm([xx(kk) - xc, yy(kk) - yc]);

    % Angle calculation matching kymograph.m methodology
    if yy(kk) > yc
        angle(kk) = acos(dot(([xx(kk) - xc, yy(kk) - yc]) / norm([xx(kk) - xc, yy(kk) - yc]), [1 0]));
    else
        angle(kk) = 2*pi - acos(dot(([xx(kk) - xc, yy(kk) - yc]) / norm([xx(kk) - xc, yy(kk) - yc]), [1 0]));
    end
end

%% --- FIRST PASS: fit r(angle) and fill middle half of RRrow ---
param = polyfit(angle, r, polyOrder);
xgrid = linspace(0, 2*pi, nBoundary);
y1 = polyval(param, xgrid);

r_theta_p   = polyder(param);
r_2theta_p  = polyder(r_theta_p);
r_theta     = polyval(r_theta_p, xgrid);
r_2theta    = polyval(r_2theta_p, xgrid);

% Curvature calculation for middle half (bins 26-75 for 101 theta bins)
for ii = nThetaBins/2 - nThetaBins/4 : nThetaBins/2 + nThetaBins/4
    sel = xgrid > theta(ii) & xgrid < theta(ii+1);
    if any(sel)
        num = ((y1(sel).^2 + r_theta(sel).^2).^(3/2));
        den = abs(y1(sel).^2 + 2*r_theta(sel).^2 - y1(sel).*r_2theta(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

RRrow(126:375) = y1(126:375);

%% --- SECOND PASS: sort and shift by pi to handle wrap-around ---
[angle2, Iord] = sort(angle);
r2 = r(Iord);

angle2 = angle2 + pi;
angle2(angle2 > 2*pi) = angle2(angle2 > 2*pi) - 2*pi;

param2 = polyfit(angle2, r2, polyOrder);
x2 = linspace(0, 2*pi, nBoundary);
y2 = polyval(param2, x2);

r_theta_p2   = polyder(param2);
r_2theta_p2  = polyder(r_theta_p2);
r_theta2     = polyval(r_theta_p2, x2);
r_2theta2    = polyval(r_2theta_p2, x2);

x2 = x2 - pi;
x2(x2 < 0) = 2*pi + x2(x2 < 0);

y2 = circshift(y2, nBoundary/2);

RRrow(1:125) = y2(1:125);
RRrow(376:500) = y2(376:500);

% Curvature for remaining quarters
for ii = 1:nThetaBins/2 - nThetaBins/4
    sel = x2 > theta(ii) & x2 < theta(ii+1);
    if any(sel)
        num = ((y2(sel).^2 + r_theta2(sel).^2).^(3/2));
        den = abs(y2(sel).^2 + 2*r_theta2(sel).^2 - y2(sel).*r_2theta2(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

for ii = nThetaBins/2 + nThetaBins/4 : nThetaBins
    sel = x2 > theta(ii) & x2 < theta(ii+1);
    if any(sel)
        num = ((y2(sel).^2 + r_theta2(sel).^2).^(3/2));
        den = abs(y2(sel).^2 + 2*r_theta2(sel).^2 - y2(sel).*r_2theta2(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

%% --- Reconstruct final smooth boundary from RRrow ---
theta_boundary = linspace(0, 2*pi, nBoundary);
XX = RRrow .* cos(theta_boundary) + xc;
YY = RRrow .* sin(theta_boundary) + yc;

% Note: caller should clamp to image bounds if needed for poly2mask
poly = [XX(:), YY(:)];

end

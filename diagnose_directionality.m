% diagnose_directionality.m
% Creates diagnostic overlay showing reference axis and tangent direction conventions
%
% Run this AFTER boundary_flows.m has processed at least one frame
% (or run the first section of boundary_flows.m up through frame 1)
%
% OUTPUT: Saves diagnostic image to outDir showing:
%   - Reference axis (theta = 0, pointing +X / right)
%   - Direction of increasing theta (CCW arrow)
%   - Sample tangent vectors at key positions (0°, 90°, 180°, 270°)
%   - Legend explaining CCW (+) vs CW (-) convention

%% ========================== CONFIGURATION ================================
% If running standalone, set these paths:
% Otherwise, these should already be defined from boundary_flows.m

if ~exist('outDir', 'var')
    error('Run boundary_flows.m first, or define outDir, visImageSeq, visPolySeq, centroidXY');
end

% Select a frame with valid data
fr_diag = find(qcFlag, 1, 'first');
if isempty(fr_diag)
    error('No frames passed QC. Run boundary_flows.m first.');
end

fprintf('Creating directionality diagnostic for frame %d...\n', fr_diag);

%% ========================== LOAD FRAME DATA ==============================
I = visImageSeq{fr_diag};
poly = visPolySeq{fr_diag};
xc = centroidXY(fr_diag, 1);
yc = centroidXY(fr_diag, 2);

% Get mean radius for visualization
R_mean = mean(RADIUS_OF_CURVATURE(fr_diag, :), 'omitnan');
if isnan(R_mean)
    R_mean = 100;  % fallback
end

%% ========================== CREATE DIAGNOSTIC FIGURE =====================
figDiag = figure('Position', [50 50 1000 900]);

% Panel 1: Full diagnostic overlay
subplot(1,1,1);
imshow(I, []); hold on;

% Plot boundary
plot(poly(:,1), poly(:,2), 'c-', 'LineWidth', 1.5);

% Plot centroid
plot(xc, yc, 'w+', 'MarkerSize', 20, 'LineWidth', 3);
plot(xc, yc, 'k+', 'MarkerSize', 18, 'LineWidth', 2);

% === REFERENCE AXIS: Theta = 0 (pointing +X, i.e., RIGHT) ===
arrow_len = R_mean * 1.3;  % extend beyond boundary
quiver(xc, yc, arrow_len, 0, 0, 'Color', 'y', 'LineWidth', 3, 'MaxHeadSize', 0.5);
text(xc + arrow_len + 20, yc, '\theta = 0°', 'Color', 'y', 'FontSize', 14, 'FontWeight', 'bold');

% === SHOW DIRECTION OF INCREASING THETA ===
% Draw curved arrow indicating direction of increasing theta
% Note: In image coords (+Y down), this appears CLOCKWISE on screen
theta_arc = linspace(0, pi/3, 50);
arc_r = R_mean * 0.5;
arc_x = xc + arc_r * cos(theta_arc);
arc_y = yc + arc_r * sin(theta_arc);  % +Y is DOWN, so this curves downward = CW visually
plot(arc_x, arc_y, 'g-', 'LineWidth', 2);
% Arrowhead at end of arc
quiver(arc_x(end-1), arc_y(end-1), arc_x(end)-arc_x(end-1), arc_y(end)-arc_y(end-1), 0, ...
    'Color', 'g', 'LineWidth', 2, 'MaxHeadSize', 2);
text(xc + arc_r*1.1, yc + arc_r*0.6, '+\theta', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');
text(xc + arc_r*1.1, yc + arc_r*0.9, '(CW on screen)', 'Color', 'g', 'FontSize', 10);

% === SAMPLE TANGENT VECTORS AT KEY ANGLES ===
% Show tangent direction at 0°, 90°, 180°, 270°
key_angles = [0, pi/2, pi, 3*pi/2];
angle_labels = {'0°', '90°', '180°', '270°'};
tangent_len = 60;  % pixels

for k = 1:numel(key_angles)
    theta_k = key_angles(k);

    % Position on boundary
    px = xc + R_mean * cos(theta_k);
    py = yc + R_mean * sin(theta_k);

    % Tangent vector (perpendicular to radial, pointing CCW)
    % CCW tangent: rotate radial by +90° → (-sin(θ), cos(θ))
    tx = -sin(theta_k);
    ty = cos(theta_k);

    % Draw tangent arrow (magenta = CCW/positive direction)
    quiver(px, py, tangent_len*tx, tangent_len*ty, 0, ...
        'Color', [1 0 1], 'LineWidth', 2, 'MaxHeadSize', 1);

    % Mark the position
    plot(px, py, 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');

    % Label
    label_offset = 25;
    text(px + label_offset*cos(theta_k), py + label_offset*sin(theta_k), ...
        angle_labels{k}, 'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.3 0.3 0.3], 'HorizontalAlignment', 'center');
end

% === ADD LEGEND / EXPLANATION BOX ===
legend_text = {
    'DIRECTIONALITY CONVENTION:'
    ''
    '  Reference: \theta = 0 points RIGHT (+X)'
    '  \theta increases CLOCKWISE on screen'
    '  (because +Y points DOWN in images)'
    ''
    '  Tangent t = (-sin\theta, cos\theta)'
    '  Points in direction of increasing \theta'
    ''
    '  v_\theta = v \cdot t  (dot product)'
    ''
    '  POSITIVE / RED: Flow CLOCKWISE (on screen)'
    '         (magenta arrows show + direction)'
    ''
    '  NEGATIVE / BLUE: Flow COUNTER-CLOCKWISE'
    '         (opposite to magenta arrows)'
};

annotation('textbox', [0.02 0.02 0.35 0.32], 'String', legend_text, ...
    'FontSize', 10, 'FontName', 'FixedWidth', ...
    'BackgroundColor', [0 0 0 0.7], 'Color', 'w', ...
    'EdgeColor', 'w', 'LineWidth', 1, ...
    'VerticalAlignment', 'top', 'Interpreter', 'tex');

% === COORDINATE SYSTEM INSET ===
% Show image coordinate system reminder
annotation('textbox', [0.72 0.02 0.26 0.18], 'String', {
    'IMAGE COORDINATES:'
    '  +X = RIGHT'
    '  +Y = DOWN (flipped!)'
    '  Origin = top-left'
    ''
    'So +theta goes CW visually'
}, 'FontSize', 9, 'FontName', 'FixedWidth', ...
    'BackgroundColor', [0.2 0.2 0.2 0.8], 'Color', 'w', ...
    'EdgeColor', 'y', 'LineWidth', 1);

title(sprintf('Directionality Diagnostic - Frame %d', fr_diag), ...
    'FontSize', 14, 'Color', 'w', 'FontWeight', 'bold');

%% ========================== SAVE DIAGNOSTIC ==============================
outFile = fullfile(outDir, 'DIAGNOSTIC_directionality_reference.png');
exportgraphics(figDiag, outFile, 'Resolution', 250);
fprintf('Saved directionality diagnostic to:\n  %s\n', outFile);

close(figDiag);

%% ========================== PRINT SUMMARY ================================
fprintf('\n');
fprintf('=======================================================\n');
fprintf('         DIRECTIONALITY REFERENCE SUMMARY\n');
fprintf('=======================================================\n');
fprintf('\n');
fprintf('COORDINATE SYSTEM:\n');
fprintf('  - Image origin: top-left corner\n');
fprintf('  - +X direction: RIGHT\n');
fprintf('  - +Y direction: DOWN (flipped from standard math!)\n');
fprintf('\n');
fprintf('ANGULAR REFERENCE (in image coords):\n');
fprintf('  - theta = 0°:   points RIGHT (+X) from centroid\n');
fprintf('  - theta = 90°:  points DOWN (+Y) from centroid\n');
fprintf('  - theta = 180°: points LEFT (-X) from centroid\n');
fprintf('  - theta = 270°: points UP (-Y) from centroid\n');
fprintf('\n');
fprintf('*** IMPORTANT: VISUAL vs MATHEMATICAL ***\n');
fprintf('  Because +Y points DOWN in image coordinates:\n');
fprintf('  - Increasing theta goes CLOCKWISE on screen\n');
fprintf('  - "CCW" in code = CLOCKWISE when viewing image\n');
fprintf('  - "CW" in code = COUNTER-CLOCKWISE when viewing image\n');
fprintf('\n');
fprintf('TANGENT VECTOR at angle theta:\n');
fprintf('  t = (-sin(theta), cos(theta))\n');
fprintf('  Points in direction of INCREASING theta\n');
fprintf('\n');
fprintf('TANGENTIAL VELOCITY v_theta:\n');
fprintf('  v_theta = dot(velocity, tangent)\n');
fprintf('  POSITIVE (+): flow in direction of increasing theta\n');
fprintf('               (CLOCKWISE when viewing image)\n');
fprintf('  NEGATIVE (-): flow opposite to increasing theta\n');
fprintf('               (COUNTER-CLOCKWISE when viewing image)\n');
fprintf('\n');
fprintf('VISUAL INTERPRETATION (looking at the image):\n');
fprintf('  RED arrows (+ velocity):  flow goes CLOCKWISE\n');
fprintf('  BLUE arrows (- velocity): flow goes COUNTER-CLOCKWISE\n');
fprintf('  Zero-crossings: convergence/divergence points (poles)\n');
fprintf('=======================================================\n');

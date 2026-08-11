% RUN_SCW_ANALYSIS  Cortex segmentation -> curvature change -> SCW strength.
%
% The pipeline, in order:
%   1. scw_cortex     segments every frame (Laplacian ridge + snap to the
%                     retardance maximum) and builds the cortical and
%                     subcortical rings around the outline.
%   2. scw_curvature  smooths each outline with piecewise polynomials in
%                     polar coordinates, computes the first principal
%                     curvature on ~2 um segments, subtracts the first
%                     frame's curvature, and quantifies SCW strength as the
%                     background-subtracted variance of the radii of
%                     curvature. Cortical/subcortical retardance is averaged
%                     in the same segments.
%   3. scw_export     overlays, movie, and per-contour CSVs.
%
% See README_SCW.md for how each step maps onto the published schema.

src = '/Users/hridaytalreja/Desktop/Mar_2026_data/2026_05_26_F9/SMS_2026_0526_1552_1/Pos0/';

R = scw_cortex(src);

% First pass without windows: look at the var(rho) panel (top right of the
% scw_curv figure) to read off when the SCW runs, then set the windows.
C = scw_curvature(R);

% Second pass with the SCW window and an equal-duration metaphase window,
% writing the CSVs next to the data:
% C = scw_curvature(R, 'ScwWindow', [30 55], 'BgWindow', [0 25], ...
%                      'OutDir', fullfile(src, 'scw_output'));

scw_export(R, src);

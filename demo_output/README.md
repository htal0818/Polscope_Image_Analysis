# Demo Output — Single-Image Analysis

Results from `run_single_image.py` on a synthetic oocyte test image.

## Input Image
![test_oocyte](test_oocyte.png)

## 1. Boundary Detection Overlay
![overlay](overlay.png)

Red = shrunk boundary | Green = Fourier-fit | Yellow = inward normals | Cyan = circle fit

## 2. Mask & Distance Transform
![mask_and_distance](mask_and_distance.png)

Left: retardance | Center: binary mask | Right: distance from cortex

## 3. Center-Out Radial Profile
![radial_profile_center_out](radial_profile_center_out.png)

## 4. Outside-In: Normal-Based
![radial_profile_normal_outside_in](radial_profile_normal_outside_in.png)

## 5. Outside-In: Distance Transform
![radial_profile_dist_outside_in](radial_profile_dist_outside_in.png)

## 6. Comparison: Normal vs Distance Transform
![comparison_normal_vs_dist](comparison_normal_vs_dist.png)

## 7. Retardance vs Angle at Cortex
![retardance_vs_angle](retardance_vs_angle.png)

## Summary

| Metric | Value |
|--------|-------|
| Oocyte radius | 151 px = 24.2 um |
| Contour retardance (mean) | 3.97 +/- 3.26 nm |
| Contour retardance (range) | 0.00 - 11.61 nm |
| Fourier harmonics | 25 |

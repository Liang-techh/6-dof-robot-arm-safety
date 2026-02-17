========================================================================
6-DOF Robotic Arm SOS Reachability Analysis & Visualization
Student: [Your Name]
Date: 2026-02-16
========================================================================

This folder contains the results of the Sum-of-Squares (SOS) based L2 
reachability analysis for a 6-DOF robotic arm.

Files:
1. work24_6dof.m       : Main analysis script (Requires YALMIP + MOSEK).
                         Computes the Lyapunov barrier function and safety sets.
2. result_6dof_arm.mat : Pre-computed results and model parameters.
3. visualize_6dof_arm.m: Interactive 3D visualization tool.
                         (Requires only standard MATLAB).

------------------------------------------------------------------------
HOW TO VIEW RESULTS (No toolboxes required):
------------------------------------------------------------------------
1. Open this folder in MATLAB.
2. Run the command: 
   >> visualize_6dof_arm

3. Interactive Controls:
   - [Start Simulation]: Plays an infinite loop of the stabilized trajectory.
   - [Inject Disturbance]: Adds random perturbations to test robustness.
   - [Stop]: Pauses simulation to allow manual slider control.
   - [Sliders]: Manually inspect the arm's safety constraints (Right Panels).

------------------------------------------------------------------------
ANALYSIS DETAILS:
------------------------------------------------------------------------
- Model: 12-state linearized dynamics (6 joints, q + dq).
- DH Parameters: Based on actual C++ kinematic definitions.
- Safety Guarantee: 
  * Eta* = 0.80 (Lyapunov Level Set)
  * Alpha* = 1.98 (Terminal Set)
  * R_init = 0.15 (Certified Initial Set Radius)

The right-side plots in the visualization show the real-time state (red dots)
remaining STRICTLY within the computed safety ellipses (black lines),
verifying the theoretical guarantees.

# 6-DOF Robotic Arm: Safe Control with SOS Verification

This project implements a **Sum-of-Squares (SOS)** based safety verification and control framework for a 6-Degree-of-Freedom (6-DOF) robotic arm. 

The core contribution is a rigorously verified safety controller that guarantees the robotic arm remains within a computed safety region (Lyapunov level set) despite disturbances and model uncertainties.

![Visualization Preview](https://github.com/user-attachments/assets/placeholder-image-url)
*(Note: You can add a screenshot here after uploading)*

## 📂 Project Structure

- **`visualize_6dof_arm.m`**: **Start Here!** Interactive 3D visualization tool.
  - Run this in MATLAB to see the arm in action.
  - Features: Infinite simulation, manual slider control, disturbance injection.
- **`result_6dof_arm.mat`**: Pre-computed safety verification results.
  - Contains the Lyapunov function parameters (`Vsym`), controller gains (`a_coef`, `c_coef`), and safety levels (`eta_star`).
- **`work24_6dof.m`**: The main analysis script (Source Code).
  - Uses **YALMIP** and **MOSEK** to solve the SOS optimization problem.
  - Generates the `.mat` results file.

## 🚀 How to Run

### Prerequisites
- **MATLAB** (Tested on R2016b+).
- *No additional toolboxes required for visualization.*
- *(Optional)* To re-run the analysis (`work24_6dof.m`), you need [YALMIP](https://yalmip.github.io/) and [MOSEK](https://www.mosek.com/).

### Quick Start
1.  Download this repository.
2.  Open MATLAB and navigate to the folder.
3.  Run the visualization command:
    ```matlab
    visualize_6dof_arm
    ```
4.  Interact with the GUI:
    - **`Start Simulation`**: Watch the certified safe trajectory.
    - **`Inject Disturbance`**: Click to test robustness.
    - **`Stop`**: Pause to manually control joints with sliders.

## 🔬 Scientific Details

- **Dynamics**: Linearized Euler-Lagrange dynamics for 6 joints (12 states).
- **Safety Verification**:
    - **Lyapunov Function**: Quadratic candidate $V(x) = x^T P x$.
    - **SOS Programming**: Used to find the largest level set $\mathcal{E} = \{x | V(x) \le \eta^*\}$ where $\dot{V}(x) < 0$ along boundaries.
    - **Result**: Certified safety radius $\eta^* = 0.80$.

## 🎥 Visualization Features

The right-side panel in the visualization tool shows real-time 2D projections of the 12D state space (e.g., $q_1$ vs $q_2$). The **red dot** (current state) is mathematically guaranteed to stay within the **black ellipse** (safety boundary).

---
*Created by [Your Name]*

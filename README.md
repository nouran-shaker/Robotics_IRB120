# ABB IRB 120 Digital Twin: Kinematics, Dynamics, and Control

## 📌 Project Overview
This repository contains the comprehensive MATLAB and Simulink codebase for modeling the 6-DOF **ABB IRB 120** industrial manipulator. Developed as the final course requirement for Industrial Robotics (MCT344/MCT342) at Ain Shams University, this project simulates an industrial pick-and-place task to cross-validate custom mathematical derivations against numerical software toolboxes.

## 🛠️ Methodologies 
The digital twin and control systems were evaluated using comparative engineering methods:
1. **Analytical Modeling:** Custom MATLAB scripts deriving the Denavit-Hartenberg (DH) kinematics, 8-branch Inverse Kinematics, Geometric Jacobian, and Euler-Lagrange dynamics from foundational mathematical principles.
2. **Robotics System Toolbox:** Numerical verification using MATLAB's built-in `rigidBodyTree` and recursive Newton-Euler inverse dynamics solvers.
3. **Simscape Multibody:** A 3D physical plant simulation utilized for trajectory validation, collision checking, and control system analysis.

## 📂 Repository Structure
* `/Analytical_Models/` - Custom closed-form MATLAB scripts (numbered for sequential execution).
* `/Toolbox_Models/` - Scripts utilizing the MATLAB Robotics System Toolbox.
* `/Comparison_Scripts/` - Automated testing scripts to calculate maximum absolute errors between the two methods.
* `/Simscape_Simulation/` - Simulink block diagrams and physical 3D models.


## 🚀 Execution Guide
Because the mathematical models and the toolbox scripts share trajectory and joint-space data, the files **must be executed in a specific order**. 

Please refer to the [RUN_INSTRUCTIONS.md](RUN_INSTRUCTIONS.md) file in this repository for the exact step-by-step guide on how to initialize the workspace and run the simulations.

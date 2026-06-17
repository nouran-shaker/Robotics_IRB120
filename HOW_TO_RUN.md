# How to Run the Project

To ensure the MATLAB Base Workspace is properly populated with the correct parameters, transformations, and joint angles, the project files **must** be executed in the exact order listed below. 

### ⚙️ Prerequisites
* **MATLAB** (R2023a or newer recommended)
* **Robotics System Toolbox**

---

### Step 1: Initialize the MATLAB Toolbox and IK Solutions
**File to run:** matlab_toolbox.m

You must start with the MATLAB Toolbox script. This file is responsible for setting up the spatial pick-and-place trajectory and generating the initial Inverse Kinematics (IK) solutions.
* **What it does:** It generates the target joint angles for *both* the numerical toolbox method and the custom analytical method, exporting these seed values to the workspace so both models follow the exact same path.

### Step 2: Run the Analytical Models (In Order)
Once the workspace is seeded with the IK solutions, navigate to the analytical scripts folder. You must run these files sequentially based on their numbering to build the mathematical model step-by-step:

1. **`1_robot_params.m`**: Initializes the DH parameters, masses, and inertias.
2. **`2_forward_kinematics.m`**: Computes the Forward Kinematics and verifies the spatial frames.
3. **`3_jacobian.m`**: Calculates the Geometric Jacobian and analyzes singularities.
4. **`4_irb120_dynamics.m`**: Computes the Inertia, Coriolis, and Gravity matrices using Euler-Lagrange equations.

### Step 3: Run the Results Comparison
**File to run:** `[compare_methods.m]`

After both the toolbox and analytical workspaces are fully populated, run the comparison script.
* **What it does:** This script aligns the two methodologies, applies the necessary parallel-axis shifts for inertia, and calculates the maximum absolute errors for the Jacobian, Inertia Matrix, Gravity Vector, and Coriolis terms. 
* **Expected Output:** The error tables will print directly to the MATLAB Command Window (expected dynamic variance is mathematically negligible, typically around `1e-11` N·m).

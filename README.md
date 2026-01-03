# Simulation of a Robust Fuzzy-PID Controller for Quadcopter Position and Attitude
In this project, a robust Fuzzy-PID controller for a 6-DOF quadcopter was designed and simulated in MATLAB & Simulink to overcome the limitations of standard PID controllers, such as instability, poor tracking, and sensitivity to uncertainties. The resulting intelligent controller demonstrated superior performance by eliminating overshoot and oscillations and significantly improving trajectory accuracy. The success of the design was quantitatively validated using MSE and ITAE performance metrics.

## 🚀 How to Run the Simulation

To run the simulation in MATLAB/Simulink, follow these steps based on the project documentation:

### 1. Prerequisites
* Install **MATLAB** (Tested on MATLAB/R2023b version).
* Ensure **Simulink** and **Control System Toolbox** are installed.

### 2. Setup
1.  Open MATLAB and navigate to the project folder (e.g., `Simulation_FuzzyPID`).
2.  **Add to Path:** Right-click on the project folder in the "Current Folder" window, select **Add to Path** > **Selected Folders and Subfolders**. This ensures all functions and blocks are accessible.

### 3. Running the Simulation
1.  **Generate Trajectories:** Run the path generation scripts (e.g., for Circle, Figure-8, or Diamond paths) to create the necessary `.mat` files (like `Path_Takeoff_Circle.mat`).
2.  **Open the Model:** Open the main Simulink file (e.g., `AC_Quadcopter_Simulation.slx`).
3.  **Load Initial Conditions:**
    * Ensure the `LOAD` block is configured to load `IC-OnGround-MotorsOff.mat` for zero initial states.
    * Or double-click the `LOAD` block to reload parameters manually.
4.  **Run:** Click the **Run** button in Simulink.

## 📊 Project Description
This project focuses on the design and implementation of a **Robust Fuzzy-PID Controller** for a 6-DOF Quadrotor. It addresses the limitations of classical PID controllers, such as poor tracking in aggressive maneuvers and sensitivity to uncertainties.

**Key Features:**
* **Simulation:** Full 6-DOF nonlinear dynamic model in Simulink.
* **Control Strategy:** Hybrid Fuzzy-PID for adaptive gain tuning.
* **Scenarios:** Tested on Circle, Figure-8, and Diamond trajectories.
* **Performance:** Validated using MSE and ITAE metrics, showing superior stability compared to standard PID.

## 📚 References & Acknowledgments
This project utilizes core dynamic models and simulation frameworks adapted from the following source:
* **Original Repository:** [dch33/quad-sim](https://github.com/dch33/quad-sim)
* **Description:** Provides the fundamental 6-DOF Quadcopter dynamic modeling and MATLAB/Simulink visualization blocks.
# LESM - Linear Elements Structure Model

LESM is a MATLAB-based structural analysis program for linear-elastic, displacement-based, static analysis of two-dimensional and three-dimensional reticulated structures. It employs the direct stiffness method and follows an Object-Oriented Programming (OOP) paradigm.

## Project Overview

- **Core Purpose:** Analysis of 2D/3D Trusses, 2D/3D Frames, and Grillages.
- **Analysis Theories:** Supports both Euler-Bernoulli (Navier) and Timoshenko beam theories for flexural behavior.
- **Implementation:** MATLAB (using OOP classes).
- **Input Format:** Neutral format files with `.lsm` extension.

### Architecture

The project is organized around several key classes:

- **`Drv` (Driver):** The central controller that manages the structural model, global stiffness matrix, and the overall analysis process.
- **`Anm` (Analysis Model):** An abstract base class for specific analysis types (e.g., `Anm_Truss2D`, `Anm_Frame3D`).
- **`Node`:** Represents structural joints and points of discretization.
- **`Elem` (Element):** An abstract base class for structural members. Subclasses like `Elem_Navier` and `Elem_Timoshenko` implement specific flexural behaviors.
- **`Lelem` (Load Element):** Manages distributed and thermal loads applied to elements.
- **`Material` & `Section`:** Define physical and geometric properties of the model.
- **`Print`:** Handles formatting and displaying analysis results.

## Getting Started

### Prerequisites

- MATLAB installed and configured on your system.

### Running the Program

The non-graphical version of LESM is executed through the `main.m` script:

1.  Open `source_code/main.m`.
2.  Specify the target model file by editing the `fileName` variable (around line 348):
    ```matlab
    fileName = 'Models/Frame2D/Frame2D_1.lsm';
    ```
3.  Run the script (F5 in MATLAB).

The results, including nodal displacements, support reactions, and element internal forces, will be printed to the MATLAB command window.

### Project Structure

- `source_code/`: Contains all `.m` source files and class definitions.
- `source_code/Models/`: Contains sample `.lsm` model files categorized by analysis type (Beam, Frame2D, Truss3D, etc.).
- `Instructions.txt`: Quick start guide for the non-graphical version.

## Development Conventions

### Coding Style
- The project uses MATLAB's `classdef` for OOP.
- Global constants are managed via `include_constants.m` and must be included in functions using them.
- Method documentation is provided within the file headers.

### Analysis Workflow
The `Drv.process()` method defines the standard analysis sequence:
1.  **Preprocessing:** Read input file and initialize objects.
2.  **DOF Setup:** Assemble the global degree-of-freedom numbering.
3.  **Assembly:** Construct the global stiffness matrix (`K`) and forcing vector (`F`).
4.  **Solution:** Solve the partitioned system of equilibrium equations ($K \cdot D = F$).
5.  **Post-processing:** Compute internal forces and displacements for each element.

### Adding New Features
- **New Analysis Model:** Inherit from the `Anm` class.
- **New Element Type:** Inherit from the `Elem` class and implement the abstract flexural methods.
- **Adding I/O:** Update `readFile.m` and `saveFile.m` to handle new tags in the `.lsm` format.

## Building and Testing
- **Building:** No compilation is required as it is a MATLAB project.
- **Testing:** Currently, testing is performed by running `main.m` against the provided sample models in the `Models/` directory. (TODO: Implement automated unit testing for core matrix operations).

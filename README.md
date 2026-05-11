# README

## Overview

This MATLAB code performs a nonlinear displacement-controlled finite element analysis of a 3D concrete portal frame using 8-node hexahedral (Hex8) elements. It calculates the global tangent stiffness matrix, applies boundary conditions, solves equilibrium iteratively using the Newton-Raphson method, computes displacement and stress/strain distributions at each load step, and generates load-displacement curves and visualizations of the deformed structure.

---

## Input

| Parameter | Description |
|-----------|-------------|
| `E0` | Initial Young's modulus of the material (30 GPa) |
| `nu` | Poisson's ratio (0.18) |
| `sigma_max` | Saturation stress of the concrete (35 MPa) |
| `mshFile` | Mesh file in Gmsh format (`GMSH/frame.msh`) |
| `fixedDOF` | Fixed support DOFs — all nodes at Z = 0 (base of columns) |
| `controlDOF` | Controlled displacement DOF — Node 318, X-direction |
| `u_target` | Target prescribed displacement (−0.058 m = −58 mm = 4% drift) |
| `nSteps` | Number of displacement increments (30) |
| `maxIter` | Maximum Newton-Raphson iterations per step (50) |
| `nrTol` | Convergence tolerance (1×10⁻³) |

---

## Output

| Output | Description |
|--------|-------------|
| `frame_output.xlsx` | Node coordinates and element connectivity exported to Excel |
| Force-Displacement Graph | Load-displacement curve plotted at each load step |
| Deformed Shape Plot | Undeformed (ghost) overlaid with deformed and stress-coloured mesh |
| Stress Contour Plots | Saturation stress, Drucker-Prager D/C ratio, principal stresses, hydrostatic pressure, and tangent modulus fields |
| Load Step Animation | 30-step animated deformation history with growing load-displacement curve |
| Console Output | Reaction forces, max/min stresses, support reactions, and displacement results |

---

## Main Steps

1. Read the Gmsh mesh file and extract node coordinates and element connectivity.
2. Initialize element and node data structures with elastic material state.
3. Apply prescribed displacement incrementally across 30 load steps.
4. At each step, iterate using Newton-Raphson until equilibrium is satisfied:
   - Compute strains and stresses at all 8 Gauss points per element
   - Update tangent modulus `E_tan` using the hyperbolic saturation material model
   - Assemble the global tangent stiffness matrix
   - Assemble the internal force vector
   - Solve for displacement correction on free DOFs only
   - Re-enforce essential boundary conditions exactly
5. Recover reaction force at the controlled DOF from the converged internal force vector.
6. Compute nodal and element-level stress and strain distributions.
7. Generate load-displacement curve, deformed shape, stress contour plots, and animation.
8. Export mesh data to Excel and print results to console.

---

## Dependencies

| File | Purpose |
|------|---------|
| `Hex8_ReadGmsh22.m` | Reads Gmsh `.msh` file (format 2.2 and 4.1) |
| `Hex8_ExtractMeshToExcel.m` | Extracts mesh data and exports to Excel |
| `Hex8_ElementStiffness.m` | Computes 24×24 element stiffness matrix |
| `Hex8_AssembleGlobalStiffness.m` | Assembles global stiffness matrix |
| `Hex8_ComputeStrainStress.m` | Computes strains, stresses and tangent modulus at Gauss points |
| `Hex8_TangentD.m` | Builds 6×6 isotropic tangent D matrix |
| `NonlinearMaterial.m` | Hyperbolic saturation material model |
| `localAssembleInternalForce.m` | Assembles internal force vector |
| `plotDeformedMesh.m` | Plots deformed shape with stress contours |
| `plotDeformedShape.m` | Plots undeformed vs deformed mesh |
| `animateLoadSteps.m` | Animates deformation across load steps |
| `MonotonicLoad.m` | Generates monotonic load step values |

---

## Material Model

Hyperbolic saturation model:

```
sigma = E0 * eps / (1 + E0*eps / sigma_max)

E_tan = E0 / (1 + E0*eps / sigma_max)^2
```

- Stress asymptotically approaches `sigma_max` with no descending branch or defined failure point
- A minimum E floor of `1e-4 × E0 = 3 MPa` is enforced to prevent the stiffness matrix from becoming singular in damaged zones

---

## Stress Fields Available for Post-Processing

| Field | Description |
|-------|-------------|
| `sigma_sat` | Saturation equivalent stress — primary failure indicator |
| `drucker_prager` | Demand/capacity ratio — D/C = 1.0 means on failure surface |
| `principal_max` | Maximum principal stress — tension cracking indicator |
| `principal_min` | Minimum principal stress — compression crushing indicator |
| `pressure` | Hydrostatic pressure p = −I₁/3 |
| `E_tan` | Tangent modulus — identifies softened and damaged zones |
| `dispMag` | Displacement magnitude |

---

## Notes

- The code assumes the input mesh file follows Gmsh format 2.2 or 4.1
- Only 8-node hexahedral (type 5) elements are read — other element types are ignored
- Tensile cracking is not explicitly modelled — tension uses the same saturation curve
- Geometric nonlinearity (P-Δ effects) is not included
- Cyclic or hysteretic loading is not supported — monotonic loading only
- The model has no defined collapse point — unlike Hognestad, Kent-Park, or Concrete Damaged Plasticity (CDP) models

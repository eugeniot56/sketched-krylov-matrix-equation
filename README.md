# Sketching-based Krylov Methods for Generalized Matrix Equations

MATLAB research code developed during the Master's thesis of Eugenio Turchet.

This repository contains Krylov-Galerkin and sketching-based projection methods for the numerical solution of generalized Lyapunov and Sylvester matrix equations, with applications to stochastic Galerkin finite element discretizations.

---

## Features

- Krylov-Galerkin solvers for generalized Lyapunov equations
- Sketching-based Krylov projection methods
- Solvers for generalized Sylvester equations
- SGFEM numerical experiments
- Comparison with MultiRB methods
- Truncated Arnoldi orthogonalization routines
- Incremental QR factorization updates
- DCT-based randomized sketching operators

---

## Repository Structure

```text
src/
    Main MATLAB solvers and utility routines

external/
    External research code used in the experiments

experiments/
    Numerical experiments and reproducibility scripts

results/
    Generated plots and numerical results

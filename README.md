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
```

The repository is organized as follows:

```text
sketched-krylov-matrix-equations/
│
├── src/
│   ├── lyapunov/
│   ├── sylvester/
│   ├── sketching/
│   └── utils/
│
├── external/
│   └── multirb/
│
├── experiments/
│   ├── lyapunov/
│   └── sgfem/
│       └── data/
│
├── results/
│
├── README.md
├── LICENSE
├── CITATION.cff
└── .gitignore
```

---

## Main Algorithms

The repository includes implementations of:

- `solve_generalized_lyapunov_krylov`
- `solve_generalized_lyapunov_sketch_krylov`
- `solve_generalized_sylvester_sketch_krylov`

along with supporting routines for:

- truncated Arnoldi orthogonalization,
- randomized sketching,
- QR factorization updates,
- reduced projected solvers.

---

## Requirements

The code was tested with:

- MATLAB R2022a or newer
- Sparse matrix support
- Signal Processing Toolbox (for DCT-based sketching)

---

## Running the Experiments

### Generalized Lyapunov Equation Experiments

Run:

```matlab
run('experiments/lyapunov/run_lyapunov_comparison.m')
```

### SGFEM Generalized Sylvester Equation Experiments

Run:

```matlab
run('experiments/sgfem/run_sgfem_comparison.m')
```

The experiment scripts automatically configure the MATLAB path using:

```matlab
addpath(genpath(...))
```

so no manual path setup is required.

---

## SGFEM Data Files

The SGFEM `.mat` files are distributed separately from the main repository.

Download the data files from the GitHub Release page and place them in:

```text
experiments/sgfem/data/
```

Expected files:

```text
TP_2.mat
TP_3.mat
TP_4.mat
TP_5.mat
TP_6.mat
```

---

## External Research Code

This repository includes and adapts external research code from:

- C.E. Powell
- D. Silvester
- V. Simoncini
- S-IFISS toolbox

Original copyrights and licenses belong to the
respective authors.

---

## References

C.E. Powell, D. Silvester, V. Simoncini,

*An Efficient Reduced Basis Solver for Stochastic Galerkin Matrix Equations*,

SIAM Journal on Scientific Computing,  
Vol. 39, No. 1, pp. A141--A163, 2017.

---

## License

This project is distributed under the BSD 3-Clause License.

See the `LICENSE` file for details.

---

## Citation

If you use this software in academic work, please cite the associated thesis and this repository.

Suggested citation:

```text
E. Turchet,
Sketching-based Krylov methods for generalized matrix equations,
Master's Thesis,
Gran Sasso Science Institute, 2026.
```

GitHub citation metadata is also provided through the `CITATION.cff` file.

---

## Reproducibility

All experiments included in the repository are intended to be fully reproducible.

The numerical scripts:
- automatically configure paths,
- load data files from the appropriate folders,
- generate convergence plots and timing comparisons.

---

## Author

Eugenio Turchet  
Gran Sasso Science Institute (GSSI)

GitHub: https://github.com/eugeniot56

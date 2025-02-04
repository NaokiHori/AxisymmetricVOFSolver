# Axisymmetric VOF Solver

[![License](https://img.shields.io/github/license/NaokiHori/AxisymmetricVOFSolver)](https://opensource.org/license/MIT)
[![Last Commit](https://img.shields.io/github/last-commit/NaokiHori/AxisymmetricVOFSolver/main)](https://github.com/NaokiHori/AxisymmetricVOFSolver/commits/main)

![cover image](https://github.com/NaokiHori/AxisymmetricVOFSolver/blob/main/cover.jpg)

A Navier-Stokes solver for interfacial two-phase flows under axisymmetric conditions.

## Dependencies

- [C Compiler](https://gcc.gnu.org)
- [GNU Make](https://www.gnu.org/software/make/)
- [OpenMP](https://www.openmp.org) (optional)

## Quick Start

Check out the readme of [the other library](https://github.com/NaokiHori/VerySimpleNSSolver).

## Methodology

We consider two incompressible and immiscible Newtonian liquids (e.g., water and air) which are separated by free and deformable surface.
The liquids are confined inside polar coordinates under axisymmetric condition.
Their dynamics are governed by the incompressibility constraint, the mass conservation, and the momentum balance.
The whole equations are solved on the general orthogonal coordinate system (see the documentation of [the other library](https://github.com/NaokiHori/SimpleTCSolver)).

The following numerical techniques are adopted to solve the equations.

- Second-order accurate central-finite-difference scheme for spatial discretization
- Predictor-corrector method (`SMAC` / fractional-step methods) for enforcing incompressibility constraint
- Volume-of-fluid (`THINC`) method for interface capturing 
- Continuum surface-tension force model for describing surface-tension force
- Approximating scheme of the variable-coefficient Poisson equation proposed by Dodd and Ferrante to facilitate the fast Poisson solver
- Energy-consistent discretization for stable integration

For more details, check-out the following pages.

- [Simple TC Solver](https://github.com/NaokiHori/SimpleTCSolver)
- [Simple Bubbly Flow Solver](https://github.com/NaokiHori/SimpleBubblyFlowSolver)


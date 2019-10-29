# Computational Physics
This repo contains 5 different projects, all written in C++. All of the codes are physics simulations.

### Overview
1. Poisson equation: this equation is solved by two different algorithms. Both of them sets up the equation as a linear equation with a tridiagonal matrix and the two solutions are found by the Thomas algorithm and by LU decomposition.
2. Buckling Beam problem: this project is essentially a quantum physics eigenvalue problem.
3. Solar system: a simulation of our own solar system, containing the Sun, all eight planets and Pluto. This is a PDE problem and was solved by using the Velocity Verlet algorithm.
4. Ising model: a thermophysics simulation of a temperature depending magnet. This project was also parallelized with MPI. The goal was to simulate various lattices and through analysis of the results extract the so-called Curie temperature.
5. Heat diffusion: the goal in this project was to study geophysics by simulating heat diffusion in the Lithosphere. The heat sources were supposedly a geological event of radioactive nuclei inserted 1300 km below surface 1 Billion years ago.

For more information regarding the projects, please see the reports directory.

#include "test_functions.h"

void test_initial_lattice() {
/*
 * Test of initial energy and magnetization of the spin_matrix initializer
 */
  int L = 2;
  int **spin_matrix = new int* [L];
  for (int spin=0; spin<L; spin++) spin_matrix[spin] = new int[L];
  double magnetization=0;
  double energy=0;
  long idum = -1;
  Initialize_spins(spin_matrix, L, true, magnetization, energy, idum);
  if (energy != -8) {
    cout << "Initial energy: " << energy << " is wrong!" << endl;
  }
  if (magnetization != 4) {
    cout << "Initial magnetization: " << magnetization << " is wrong!" << endl;
  }
}

void test_energy_diff() {
/*
 * Test of the energy difference function in metropolis.cpp
 */
  int L = 3;
  int **spin_matrix = new int* [L];
  for (int spin=0; spin<L; spin++) spin_matrix[spin] = new int[L];
  for (int i=0; i<L; i++) {
    for (int j=0; j<L; j++) {
      spin_matrix[i][j] = 1;
    }
  }
  spin_matrix[1][0] = -1;
  int d_energy = energy_diff(1,1,L,spin_matrix);
  if (d_energy != 4) {
    cout << "energy_diff() function in metropolis.cpp is wrong!" << endl;
  }
}

#include "spin_initializer.h"

void Initialize_spins(int** spin_matrix, int L, bool order, double& magnetization,
  double& energy, long& idum) {
/*
 * Initialize a Lattice of size LxL, with either ordered configuration
 * (all spins pointing up) or a random configuration. It also finds the
 * initial energy and magnetization
 */
  if (order==true) { // ordered matrix elements equal to 1
    magnetization = L*L;
    for (int xs=0; xs<L; xs++) {
      for (int ys=0; ys<L; ys++) {
        spin_matrix[xs][ys] = 1;
      }
    }
  } else { // initialize a random-spin matrix
    for (int xs=0; xs<L; xs++) {
      for (int ys=0; ys<L; ys++) {
        int r = (int) (ran2(&idum)*L)%2; // either 0 or 1
        //(int) (ran1(&idum)*L)
        if (r==0) {
          spin_matrix[xs][ys] = -1;
          magnetization--;
        }
        else {
          spin_matrix[xs][ys] = 1;
          magnetization++;
        }
      }
    }
  }
  // calculate the initial energy
  for (int ix=0; ix<L; ix++) {
    for (int iy=0; iy<L; iy++) {
      // std::cout << spin_matrix[ix][iy] << " ";
      int mid   = spin_matrix[ix][iy];
      int down  = spin_matrix[ix][(iy+L+1)%(L)];
      int right = spin_matrix[(ix+L-1)%(L)][iy];
      energy -= mid*(down + right);
    }
    // std::cout << "\n";
  }
  // std::cout << "initial energy: " << energy << "\n";
}

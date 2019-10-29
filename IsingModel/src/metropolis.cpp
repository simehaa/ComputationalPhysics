#include "metropolis.h"

void metropolis(int** spin_matrix, int L, double& energy, double& magnetization,
                int& acceptedConfigs, double w[17], long& idum) {
  /*
   * Performs on Monte Carlo sweep of the lattice.
   * Random indices [ix][iy] are chosen each iteration
   * Energy difference (when spin[ix][iy] is flipped) is found using energy_diff()
   * Accept or decline move according to Metropolis algorithm
   * 1 Call to metropolis calls ran1 3xLxL times
   */
  int N = L*L;
  for (int i=0; i<N; i++) {
    // choose random indices
    int ix = (int) (ran2(&idum)*L);
    int iy = (int) (ran2(&idum)*L);
    int d_energy = energy_diff(ix,iy,L,spin_matrix);
    if (w[d_energy+8] >= ran2(&idum)) {
      energy += (double) d_energy;
      acceptedConfigs++;
      spin_matrix[ix][iy] *= -1;
      magnetization += (double) spin_matrix[ix][iy]*2;
    }
  }
}

int energy_diff(int ix, int iy, int L, int** spin_matrix) {
  /*
  Calulate energy difference if spin at position ix,iy is flipped.
  Takes care of periodic boundary conditions
  */
  int mid   = spin_matrix[ix][iy];
  int up    = spin_matrix[ix][(iy+L-1)%(L)];
  int down  = spin_matrix[ix][(iy+L+1)%(L)];
  int left  = spin_matrix[(ix+L+1)%(L)][iy];
  int right = spin_matrix[(ix+L-1)%(L)][iy];
  return 2*mid*(up + down + left + right);
}

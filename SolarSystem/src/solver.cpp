#include "solver.h"

solver::solver(){}
solver::solver(double DT, double TotalTime, std::vector<planet> Planets_List, bool sunfixed)
{
  dt = DT;
  totaltime = TotalTime;
  planets_list = Planets_List;
  number_of_planets = planets_list.size();
  N = totaltime/dt;
  if (sunfixed) start = 1; else start = 0;
}

void solver::velocity_verlet_solve(cube& positional_tensor)
{
  double hh_2 = dt*dt/2.0;
  double h_2 = dt/2.0;
  // save initial positions in a positional tensor:
  for (int j=start; j<number_of_planets; j++){
    positional_tensor(0, 0, j) = planets_list[j].x;
    positional_tensor(1, 0, j) = planets_list[j].y;
    positional_tensor(2, 0, j) = planets_list[j].z;
  }
  // System's initial energy
  /*double etot_i;
  for (int j=0; j<number_of_planets; j++){
    etot_i += planets_list[j].kinetic_energy();
    for (int i=0; i<number_of_planets; i++) {
      if (j != i) {
        double r2 = planets_list[j].distance(planets_list[i]);
        etot_i += planets_list[j].potential_energy(r2, planets_list[i]);
      }
    }
  }
  cout << "Systems's initial energy: " << etot_i << endl;*/

  // Integration loop
  for (int t=0; t<N-1; t++) {
    /* Matrices acceleration_matrix_old, acceleration_matrix_new
    fill store acceleration components (ax, ay, az) for all planets which
    is the sum of all contributions from all other planets */
    mat acceleration_matrix_old = zeros<mat>(3,number_of_planets); // store a_t-1
    mat acceleration_matrix_new = zeros<mat>(3,number_of_planets); // store a_t
    find_acc_for_all_planets(acceleration_matrix_old); // fill a_t-1 matrix

    for (int j=start; j<number_of_planets; j++) // start j at 1 for Sun fixed
    {
      // update x with a second degree Taylor polynomial
      positional_tensor(0, t+1, j) = planets_list[j].x + dt*planets_list[j].vx + hh_2*acceleration_matrix_old(0, j);
      positional_tensor(1, t+1, j) = planets_list[j].y + dt*planets_list[j].vy + hh_2*acceleration_matrix_old(1, j);
      positional_tensor(2, t+1, j) = planets_list[j].z + dt*planets_list[j].vz + hh_2*acceleration_matrix_old(2, j);
    }

    for (int j=start; j<number_of_planets; j++)
    {
      // update positions
      planets_list[j].x = positional_tensor(0, t+1, j);
      planets_list[j].y = positional_tensor(1, t+1, j);
      planets_list[j].z = positional_tensor(2, t+1, j);
    }

    // find acceleration components at step t
    find_acc_for_all_planets(acceleration_matrix_new);

    for (int j=start; j<number_of_planets; j++)
    {
      // update v
      planets_list[j].vx += h_2*(acceleration_matrix_new(0, j) + acceleration_matrix_old(0, j));
      planets_list[j].vy += h_2*(acceleration_matrix_new(1, j) + acceleration_matrix_old(1, j));
      planets_list[j].vz += h_2*(acceleration_matrix_new(2, j) + acceleration_matrix_old(2, j));
    }
  }
  // Systems final energy
  /*double etot_f;
  for (int j=0; j<number_of_planets; j++){
    etot_f += planets_list[j].kinetic_energy();
    for (int i=0; i<number_of_planets; i++) {
      if (j != i) {
        double r2 = planets_list[j].distance(planets_list[i]);
        etot_f += planets_list[j].potential_energy(r2, planets_list[i]);
      }
    }
  }
  cout << "Systems's final energy: " << etot_f << endl;*/
}

void solver::find_acc_for_all_planets(mat& acceleration_matrix) {
  // p is index for current planet and op is index for other planet
  // find acceleration contribution from all other planets
  double r;
  vec acc;
  for (int p=0; p<number_of_planets; p++) {
    for (int op=0; op<number_of_planets; op++) {
      // avoid finding acceleration contribution from current planet
      if (op != p) {
        r = planets_list[p].distance(planets_list[op]);
        acc = planets_list[p].acceleration(r, planets_list[op]); // acceleration vector
        acceleration_matrix(0, p) += acc[0];
        acceleration_matrix(1, p) += acc[1];
        acceleration_matrix(2, p) += acc[2];
      }
    }
  }
}

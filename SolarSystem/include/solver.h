#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>
#include <vector>
#include <cmath>
#include "planet.h"

using namespace arma;
using namespace std;

class solver
{
public:
  // initial values
  double dt;
  double totaltime;
  int N;
  int number_of_planets;
  int start; // if first element in planet_list is desired fixed at origin start=1, else 0
  std::vector<planet> planets_list;
  cube positional_tensor;

  // constructor
  solver();
  solver(double DT, double TotalTime, std::vector<planet> Planets_List, bool sunfixed);

  // functions
  void velocity_verlet_solve(cube& positional_tensor);
  void find_acc_for_all_planets(mat& acceleration_matrix);
};

#endif // SOLVER_H

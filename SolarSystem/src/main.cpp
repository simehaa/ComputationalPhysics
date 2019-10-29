#include <iostream>
#include <vector>
#include <armadillo>
#include "initializer.h"
#include "solver.h"
#include "planet.h"
#include "write_file.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  /*
  Units used:
  1 unit of time = 1 yr
  1 unit of length = 1 AU
  1 unit of velocity = 1 AU/yr
  1 unit of mass = 1.99*10^30 kg = 1 mass_sun
  */
  double T = 248; // one Pluto year
  double dt = 0.001; // Should not be larger than 0.001
  int N = T/dt;
  vec time_vec = linspace<vec>(0,T,N);

  // read input variables and create a vector containing planet objects
  vector<planet> planets_list;
  planets_list = init_planet_list("./data/initial_values.txt");
  int number_of_planets = planets_list.size();
  cout << "Number of integration points = " << N << endl;
  cout << "Number of planets = " << number_of_planets << endl;

  // initialize solver instance
  solver VelVerlet_Solarsystem;
  bool sunfixed = false;
  VelVerlet_Solarsystem = solver(dt, T, planets_list, sunfixed);
  cube positional_tensor = zeros<cube>(3,N,number_of_planets);
  cout << "Simulating Solar System..." << endl;
  // solve: most time consuming part of program
  VelVerlet_Solarsystem.velocity_verlet_solve(positional_tensor);
  // write positions to files
  write_new_file(positional_tensor, time_vec, planets_list); // set to only write every 1000th data point
  return 0;
}

/* Run example
simen@simen-ubuntu:~/steinngithub/FYS3150/Project3/OO_solar_system$ ./test.x
Number of integration points = 248000
Number of planets = 10
Simulating Solar System...
File "./data/planet0.txt" written.
File "./data/planet1.txt" written.
File "./data/planet2.txt" written.
File "./data/planet3.txt" written.
File "./data/planet4.txt" written.
File "./data/planet5.txt" written.
File "./data/planet6.txt" written.
File "./data/planet7.txt" written.
File "./data/planet8.txt" written.
File "./data/planet9.txt" written.
*/
/* Three body run:
simen@simen-ubuntu:~/steinngithub/FYS3150/Project3/OO_solar_system$ ./test.x
Number of integration points = 15000
Number of planets = 3
Earth's initial energy: -39.4302
Jupiter's initial energy: -7.33643
Simulating Solar System...
Earth's final energy: -39.4302
Jupiter's final energy: -7.33643
*/

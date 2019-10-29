#ifndef PLANET_H
#define PLANET_H
#include <armadillo>
#include <math.h> // for M_PI
#include <string>
#include <cmath> // for sqrt 

using namespace std;
using namespace arma;

class planet
{
public:
  // initial values and constants
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  double mass;
  string name;
  double gravity_factor = 4*M_PI*M_PI;
  double c = 63239.7263; // au/year
  double cc = c*c;

  // constructor
  planet();
  planet(double Mass, double X, double Y, double Z, double VX, double VY, double VZ);

  // functions
  double distance(const planet& otherplanet);
  vec acceleration(double r, const planet& otherplanet); // acceleration in x,y,z direction
  double potential_energy(double r, const planet& otherplanet); //Only calculating the potential in relation to the Sun.
  double kinetic_energy();
  void setName(string Name);
};

#endif // PLANET_H

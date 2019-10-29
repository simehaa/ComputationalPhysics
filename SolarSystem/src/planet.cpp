#include "planet.h"

planet::planet(){}
planet::planet(double Mass, double X, double Y, double Z, double VX, double VY, double VZ) {
  mass = Mass;
  x = X;
  y = Y;
  z = Z;
  vx = VX;
  vy = VY;
  vz = VZ;
}

double planet::distance(const planet& otherplanet) {
  double xdiff = x - otherplanet.x;
  double ydiff = y - otherplanet.y;
  double zdiff = z - otherplanet.z;
  return sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);
}

vec planet::acceleration(double r, const planet& otherplanet) {
  double GR_factor = 1;
  // GENERAL RELATIVITY extension
  // cross product squared:
  /*
  double c1 = y*vz - z*vy;
  double c2 = z*vx - x*vz;
  double c3 = x*vy - y*vx;
  double ll = c1*c1 + c2*c2 + c3*c3;
  GR_factor += 3*ll/(cc*r*r);
  */
  // End of general relativity extension
  double G = - GR_factor * gravity_factor * otherplanet.mass / (r*r*r);
  vec acc = zeros<vec>(3);
  // treat other planet as origin, to find x, y and z component: i.e. (x/r * G)
  acc[0] = (x - otherplanet.x)*G;
  acc[1] = (y - otherplanet.y)*G;
  acc[2] = (z - otherplanet.z)*G;
  return acc;
}

double planet::potential_energy(double r, const planet& otherplanet) {
  return - 4*M_PI*M_PI*mass*otherplanet.mass/r;
}

double planet::kinetic_energy() {
  return 0.5*mass*(vx*vx + vy*vy + vz*vz);
}

void planet::setName(string Name) {
  name = Name;
}

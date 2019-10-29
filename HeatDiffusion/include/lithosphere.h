#ifndef LITHOSPHERE_H
#define LITHOSPHERE_H
#include <armadillo>
#include <string>
#include <fstream>
#include <math.h>

using namespace arma;
using namespace std;

void JacobiMethodLithosphere(int situation);
void boundaryNoHeat(mat& u, int nx, int ny);
void boundaryNaturalHeat(mat& u, int nx, int ny);

#endif // LITHOSPHERE_H

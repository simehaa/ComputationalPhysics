#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "ran2.h"
using namespace std;

void metropolis(int** spin_matrix, int L, double& energy, double& magnetization,
                int& acceptedConfigs, double w[17], long& idum);
int energy_diff(int ix, int iy, int L, int** spin_matrix);

#endif // METROPOLIS_H

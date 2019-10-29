#ifndef SPIN_INITIALIZER_H
#define SPIN_INITIALIZER_H
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "ran2.h"

void Initialize_spins(int** spin_matrix, int L, bool order, double& magnetization, double& energy, long& idum);

#endif // SPIN_INITIALIZER_H

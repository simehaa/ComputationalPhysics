#ifndef WRITE_FILE_H
#define WRITE_FILE_H
#include <armadillo>
#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>
#include "planet.h"

using namespace std;
using namespace arma;

void write_new_file(cube positions, vec time, vector<planet> planets_list);

#endif // WRITE_FILE_H

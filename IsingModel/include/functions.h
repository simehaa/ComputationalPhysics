#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <time.h>
#include <cmath>
#include "spin_initializer.h"
#include "metropolis.h"

using namespace std;

void lattice_solve_2x2(int max_MC_steps, double temp, long idum);
void lattice_calculator(int L, int MC_steps, long idum);
void accepted_configs(int L, int MC_steps, long idum);
vector<int> prob_distribution(vector<int> energy_vec);
void phase_transition(int L, double temp, int equiltime, double &e_avg,
                      double &e2_avg, double &m_avg, double &m2_avg,
                      int MC_steps, long& idum);
void write_bin_file_double(string fn, vector<double> write_vec);
void write_bin_file_int(string fn, vector<int> write_vec);
double* w_array(double temp);
int equilibrium_time(int L, double temp, long &idum, double* w);
void equilibrium_time_distribution(int L, double temp, long &idum, int samples);

#endif // FUNCTIONS_H

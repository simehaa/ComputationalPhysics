#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <mpi.h>
#include <time.h>

#include "metropolis.h"
#include "spin_initializer.h"
#include "test_functions.h"
#include "functions.h"

using namespace std;

int main(int argc, char* argv[]) {
  long idum = -11; // RNG seed requires a negative integer

  // finds expectation values for 10^1, 10^3, 10^5 ,... MC cycles:
  // lattice_solve_2x2(10000000,1.0,idum);

  // test functions:
  // test_initial_lattice();
  // test_energy_diff();

  // calculates the mean values and probability distribution
  // lattice_calculator(20, 10000, idum);

  // calculates the equilibrium time for a given spin matrix, temperature & seed
  // equilibrium_time_distribution(20, 1.0, idum, 1000000);

  // calculates the average of total number of accepted configurations as a
  // function of Mc cycles.. Also calculates the total number of accepted
  // configurations per MC cycle per spin as a function of the temperature.
  // accepted_configs(20, 100000, idum);

  // phase transition, parallelized using MPI:
  // finds <E>, <M>, C_V and chi for different L and T
  /*
  int numprocs, my_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int MCS = 1000000;
  ofstream ofile;
  bool fileBool=true;
  double start, loopstart;
  start = clock();
  // Lattice loop
  for (int L=40; L<=40; L+=20)
  {
    if (my_rank==0) cout << "-----\nL: " << L << "\n";
    int MC_steps = MCS/numprocs;
    int equiltime = 250*L;
    long idum = - 10 - my_rank;
    double data_vec[4]={0};
    double allocate[4]={0};

    vector<double> temp_vec;
    //for (double t=2.1; t<=2.2; t+=0.05) temp_vec.push_back(t);
    //for (double t=2.21; t<=2.25; t+=0.01) temp_vec.push_back(t);
    //for (double t=2.26; t<=2.28; t+=0.005) temp_vec.push_back(t);
    //for (double t=2.29; t<=2.32; t+=0.01) temp_vec.push_back(t);
    //for (double t=2.37; t<=2.42; t+=0.05) temp_vec.push_back(t);
    temp_vec.push_back(2.50); temp_vec.push_back(2.60);

    // initialize spin_matrix for lattice size
    int **spin_matrix = new int* [L];
    for (int spin=0; spin<L; spin++) spin_matrix[spin] = new int[L];
    double magnetization=0, energy=0;
    Initialize_spins(spin_matrix, L, false, magnetization, energy, idum);
    double temp;
    int N = temp_vec.size();

    // Temperature loop
    for (int i=0; i<N; i++)
    {
      loopstart = clock();
      // Initialization
      temp = temp_vec[i];
      double *w;
      w = w_array(temp);
      int acceptedConfigs=0;
      // Reach equilibrium for first temperature
      if (i == 0)
      {
        for (int mc=0; mc<equiltime; mc++)
        {
          metropolis(spin_matrix,L,energy,magnetization,acceptedConfigs,w,idum);
        }
      }
      // MC cycles
      double e_avg=0, e2_avg=0, m_avg=0, m2_avg=0;
      for (int mc=0; mc<MC_steps; mc++)
      {
        e_avg += energy;
        e2_avg += energy*energy;
        m_avg += fabs(magnetization);
        m2_avg += magnetization*magnetization;
        metropolis(spin_matrix,L,energy,magnetization,acceptedConfigs,w,idum);
      }

      data_vec[0] = e_avg/MC_steps;
      data_vec[1] = e2_avg/MC_steps;
      data_vec[2] = m_avg/MC_steps;
      data_vec[3] = m2_avg/MC_steps;

      for (int i=0; i<4; i++)
      {
        MPI_Reduce(&data_vec[i], &allocate[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      }

      if (my_rank==0)
      {
        cout << "Temperature T=" << temp;
        for (int i=0; i<4; i++)
        {
          allocate[i] /= numprocs; //Normalize properly.
          allocate[i] /= L*L; //Make it "Per atom spin"
        }
      }

      //Write files traditionally. Only open the file for the first temperature
      //and close the file once the temperature loop is completed. This is done
      //for every lattice size desired.

      /*
      if (my_rank==0 && fileBool==true)
      {
        string filename = "./data/lattice_";
        filename.append(to_string(L) + ".txt");
        ofile.open(filename, std::ofstream::out | std::ofstream::trunc);
        ofile << setw(20) << "T:" << setw(20) << "E:" << setw(20) << "E^2: ";
        ofile << setw(20) << "M:" << setw(20) << "M^2: " << endl;

        ofile << setw(20) << setprecision(10) << temp;
        ofile << setw(20) << setprecision(10) << allocate[0];
        ofile << setw(20) << setprecision(10) << allocate[1];
        ofile << setw(20) << setprecision(10) << allocate[2];
        ofile << setw(20) << setprecision(10) << allocate[3] << endl;

        fileBool=false;
      }
      else if (my_rank==0)
      {
        ofile << setw(20) << setprecision(10) << temp;
        ofile << setw(20) << setprecision(10) << allocate[0];
        ofile << setw(20) << setprecision(10) << allocate[1];
        ofile << setw(20) << setprecision(10) << allocate[2];
        ofile << setw(20) << setprecision(10) << allocate[3] << endl;
      }
      if (my_rank==0)
      {
        double timeElapsed = (clock()-loopstart)/CLOCKS_PER_SEC;
        cout << ", time: " << timeElapsed << " s." << endl;
      }
    } // temperature loop end
    //ofile.close();
    fileBool=true;
    if (my_rank==0)
    {
      double timeElapsed = (clock()-start)/CLOCKS_PER_SEC;
      cout << "L: " << L << " calculation completed after " << timeElapsed << "s." << endl;
    }
    for(int i=0; i<L; ++i) delete[] spin_matrix[i]; delete[] spin_matrix;
  } // L loop end

  MPI_Finalize();*/
  return 0;
}

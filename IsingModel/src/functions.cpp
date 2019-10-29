#include "functions.h"

double* w_array(double temp) {
 /*
  * Function to produce an array containing the probability ratios for the
  * five different energy differences. The array is indexed from 0 to 17,
  * But the probability ratios is at index 0, 4, 8, 12, 16 and represent
  * energy difference -8, -4, 0, +4, +8
  */
  double beta = 1./temp;
  static double w[17];
  for (int i=0; i<17; i++) {
    if (i%4==0) {
      w[i] = exp(-beta*(i-8));
    } else w[i] = 0;
  }
  return w;
}

int equilibrium_time(int L, double temp, long &idum, double *w){
/*
 * Does a simulation of an LxL lattice until energy is at the ground state,
 * returns the number of MC cycles required for equilibration
 */

  int **spin_matrix = new int* [L];
  for (int spin=0; spin<L; spin++) spin_matrix[spin] = new int[L];
  double *w;
  w = w_array(temp);
  double magnetization=0, energy=0;
  int acceptedConfigs=0;
  Initialize_spins(spin_matrix, L, false, magnetization, energy, idum);

  int S=0, maxtime=3000; // upper boundary to the simulation.
  int spins = L*L;

  while (S<=maxtime){
    S++;
    metropolis(spin_matrix,L,energy,magnetization,acceptedConfigs,w,idum);
    if (energy/spins==-2){
      break;
    }
  }
  for(int i=0; i<2; ++i) delete[] spin_matrix[i]; delete[] spin_matrix;
  return S;
}

void equilibrium_time_distribution(int L, double temp, long &idum, int samples) {
/*
 * Calls the function equilibrium_time() <samples> number of times.
 * Counts how many MC cycles required for each run to reach equilibrium
 * Creates a histogram and writes binary file with data
 */
  vector<int> S_histogram; vector<int> S_vec;
  double mean=0;
  int S, maxEquil=3300; // maxEquil represents the largest value of S.
  // fill the S vector with equilibrium times.

  double *w;
  w = w_array(temp);

  double start, finish;
  start = clock();
  for (int i=0; i<samples; i++){
    S = equilibrium_time(20, 1.0, idum, w);
    S_vec.push_back(S);
    mean+=S;
  }
  mean /= samples;
  cout << "Mean: " << mean << endl;

  // loop over all possible equilibrium times
  for (int i=0; i<maxEquil; i++){
    int counter=0;
    for (int j=0; j<samples; j++) if (S_vec[j]==i) counter++;
    S_histogram.push_back(counter);
  }
  finish = clock();
  double timeElapsed = (finish-start)/CLOCKS_PER_SEC;
  cout << "Time taken " << timeElapsed << "s\n\n";

  string fn;
  fn = "S_distribution.bin";
  write_bin_file_int(fn, S_histogram);
}

void write_bin_file_double(string fn, vector<double> write_vec){
  /*
   * Function which writes vectors of type <double> as binary files.
   */
  ofstream ofile;
  ofile.open("data/" + fn, ofstream::binary);
  ofile.write(reinterpret_cast<const char*> (write_vec.data()),write_vec.size()*sizeof(double));
  ofile.close();
  cout << "File " << fn << " written." << "\n";
}

void write_bin_file_int(string fn, vector<int> write_vec){
  /*
   * Function which writes vectors of type <int> as binary files
   */
  ofstream ofile;
  ofile.open("data/" + fn, ofstream::binary);
  ofile.write(reinterpret_cast<const char*> (write_vec.data()),write_vec.size()*sizeof(int));
  ofile.close();
  cout << "File " << fn << " written." << "\n";
}

void lattice_solve_2x2(int max_MC_steps, double temp, long idum) {
  /*
   * Calculate expectation values for Energy and Magnetism, Heat capacity
   * and magnetic susceptibility for an LxL lattice for a large number of
   * Monte Carlo cycles and temperature equal to T=1.0. The average of ten
   * simulations with different seeds is used to avoid random biases.
   */

  // Outer loop with various number of MC cycles: 1e1, 1e3, 1e5, 1e7
  for (int MC_steps=10; MC_steps<=max_MC_steps; MC_steps*=100) {
    // Start time and create avg10arr which stores data
    double start, finish;
    start = clock();
    double avg10arr[10][6]; // array to calculate the average of ten runs.
    int **spin_matrix = new int* [2];
    for (int spin=0; spin<2; spin++) spin_matrix[spin] = new int[2];

    // loop over 10 sample runs to generate average values
    for (int avg10=0; avg10<10; avg10++) {
      // initialization
      double magnetization=0, energy=0;
      Initialize_spins(spin_matrix, 2, false, magnetization, energy, idum);
      double *w;
      w = w_array(temp);
      vector<int> energy_vec, magnet_vec;
      int acceptedConfigs=0;
      double e_avg=0,e2_avg=0,m_avg=0,m2_avg=0,cv=0,chi=0;

      // Monte Carlo cycles
      for (int mc=0; mc<MC_steps; mc++) {
        energy_vec.push_back(energy);
        magnet_vec.push_back(fabs(magnetization));
        e_avg += energy;
        e2_avg += energy*energy;
        m_avg += fabs(magnetization);
        m2_avg += magnetization*magnetization;
        acceptedConfigs = 0;
        metropolis(spin_matrix,2,energy,magnetization,acceptedConfigs,w,idum);
      }

      // Save averages over 10 runs
      avg10arr[avg10][0] = (double) e_avg/MC_steps;
      avg10arr[avg10][1] = (double) e2_avg/MC_steps;
      avg10arr[avg10][2] = (double) m_avg/MC_steps;
      avg10arr[avg10][3] = (double) m2_avg/MC_steps;
      avg10arr[avg10][4] = avg10arr[avg10][1] - avg10arr[avg10][0]*avg10arr[avg10][0]; // <E^2> - <E>^2
      avg10arr[avg10][5] = avg10arr[avg10][3] - avg10arr[avg10][2]*avg10arr[avg10][2]; // <M^2> - <M>^2
    }
    for(int i=0; i<2; ++i) delete[] spin_matrix[i]; delete[] spin_matrix;

    // Compute averages
    double avgs[6];
    for (int var=0; var<6; var++){
      double val=0;
      for (int run=0; run<10; run++){
        val+=avg10arr[run][var];
      }
      avgs[var] = val/10;
    }

    // Analytic values to compute relative error
    double analytics[6] = {-7.9839, 63.8714, 3.9946, 15.9732, 0.1283, 0.016};
    cout << "10 Simulations for " << MC_steps << " MC cycles produced expectation values:" << endl;
    cout << setprecision(8) << "E:   " << avgs[0] << ". Error= " << fabs(avgs[0]-analytics[0])/analytics[0] << endl;
    cout << setprecision(8) << "E^2: " << avgs[1] << ". Error= " << fabs(avgs[1]-analytics[1])/analytics[1] << endl;
    cout << setprecision(8) << "M:   " << avgs[2] << ". Error= " << fabs(avgs[2]-analytics[2])/analytics[2] << endl;
    cout << setprecision(8) << "M^2: " << avgs[3] << ". Error= " << fabs(avgs[3]-analytics[3])/analytics[3] << endl;
    cout << setprecision(8) << "Cv:  " << avgs[4] << ". Error= " << fabs(avgs[4]-analytics[4])/analytics[4] << endl;
    cout << setprecision(8) << "chi: " << avgs[5] << ". Error= " << fabs(avgs[5]-analytics[5])/analytics[5] << endl;

    finish = clock();
    double timeElapsed = (finish-start)/CLOCKS_PER_SEC;
    cout << "Time taken " << timeElapsed << "\n\n";
  }
}

void lattice_calculator(int L, int MC_steps, long idum) {
  /*
   * Analyze the evolution of a 20x20 lattice for two temperatures
   * T=1.0 and T=2.4 and two initial states ordered/random.
   * Write files with: MC-cycles, Energy and Magnetization.
   */

  int **spin_matrix = new int* [L];
  for (int spin=0; spin<L; spin++) spin_matrix[spin] = new int[L];
  // loop over ordered and random initial state
  for (int order=0; order<2; order++)
  {
    // loop over two temperatures T=1.0 and T=2.4
    for (double temp=1.0; temp<=2.4; temp+=1.4)
    {
      // initialization
      double magnetization=0, energy=0;
      int acceptedConfigs=0;
      vector<int> energy_vec, magnet_vec, mc_cycles_vec;
      double *w;
      w = w_array(temp);
      Initialize_spins(spin_matrix, L, order, magnetization, energy, idum);

      // Monte Carlo cycles
      for (int mc=0; mc<MC_steps; mc++) {
        mc_cycles_vec.push_back(mc);
        energy_vec.push_back(energy);
        magnet_vec.push_back(magnetization);
        /*if (energy/400. == -2.0) {
          cout << "Equilibrium reached at mc step: " << mc << endl; break;
        }*/
        metropolis(spin_matrix,L,energy,magnetization,acceptedConfigs,w,idum);
      }

      // calculate the probability distribution, mean and std of the results
      vector<int> probvec = prob_distribution(energy_vec);

      // Write the data to binary files accordingly
      string fn_E = "equiltime_"+to_string(order)+"_T"+to_string((int)temp)+ "_E.bin";
      string fn_M = "equiltime_"+to_string(order)+"_T"+to_string((int)temp)+ "_M.bin";
      string fn_P = "Eprob_"+to_string(order)+"_T"+to_string((int)temp)+".bin";
      string fn_MC = "equiltime_MC.bin";
      write_bin_file_int(fn_E,energy_vec);
      write_bin_file_int(fn_M,magnet_vec);
      write_bin_file_int(fn_MC,mc_cycles_vec);
      // write_bin_file_int(fn_P,probvec);
    }
  }
  for(int i=0; i<L; ++i) delete[] spin_matrix[i]; delete[] spin_matrix;
}

void accepted_configs(int L, int MC_steps, long idum) {
  /*
   * Calculate number of accepted configurations as a function of temperature
   * and as a function of MC cycles
   */
  int **spin_matrix = new int* [L];
  for (int spin=0; spin<L; spin++) spin_matrix[spin] = new int[L];
  vector<double> temp_vec, accepted_vec_T;
  vector<int> mc_vec, accepted_vec_MC;

  // loop over temperatures
  for (double temp = 1.0; temp<=3.0; temp+=0.05)
  {
    // initialization
    double magnetization, energy;
    double *w;
    w = w_array(temp);
    Initialize_spins(spin_matrix, L, false, magnetization, energy, idum);

    // Monte Carlo cycles
    int acceptedConfigs=0;
    for (int mc=0; mc<MC_steps; mc++)
    {
      metropolis(spin_matrix,L,energy,magnetization,acceptedConfigs,w,idum);
      if (temp==1.0)
      {
        accepted_vec_MC.push_back(acceptedConfigs);
        mc_vec.push_back(mc);
      }
    }

    // average accepted configs as a function of T. This should not be sampled
    // until equilibrium is reached. Loop first through equiltime before
    // sampling the MC per T.

    double acc_T_perMC = (double)acceptedConfigs/MC_steps/L/L;
    accepted_vec_T.push_back(acc_T_perMC);
    temp_vec.push_back(temp);
  }
  for(int i=0; i<L; ++i) delete[] spin_matrix[i]; delete[] spin_matrix;

  // write the arrays of type <double> to binary files.
  vector<double> vec_list1[2] = {temp_vec, accepted_vec_T};
  string string_list1[2] = {"acc_temps.bin", "acceptedconfigs_T.bin"};
  for (int i=0; i<2; i++){
    write_bin_file_double(string_list1[i], vec_list1[i]);
  }

  // write the arrays of type <int> to binary files.
  vector<int> vec_list2[2] = {mc_vec, accepted_vec_MC};
  string string_list2[2] = {"acc_mc.bin", "acceptedconfigs_MC.bin"};
  for (int i=0; i<2; i++){
    write_bin_file_int(string_list2[i], vec_list2[i]);
  }
}

vector<int> prob_distribution(vector<int> energy_vec) {
 /*
  * Function which computes the probability function P(E) of a lattice L=20,
  * for both temperatures T=1.0 and T=2.4 and initial state random/ordered.
  * The probability is found by counting the number of times a given energy
  * appears in the computations from previous results. This is done after the
  * steady state is reached. These results are then compared with the computed
  * variance in energy.
  */

  int total = energy_vec.size();
  // decide where to start counting. Should be after equilibrium. The seed
  // idum=-6 causes an equilibrium time of 2300 MC cycles. Total*0.03 should
  // account for avoiding the first three thousand MC cycles.
  int startVal = total*0.03;
  double mean, std, sum=0, variance=0, count=0, newSize;
  newSize = total-startVal;

  // calculate the mean and standard deviation for energy analyses.
  for (int i=startVal; i<total; i++) sum+=energy_vec[i];
  mean = sum/newSize;
  for (int i=startVal; i<total; i++) variance+=pow(energy_vec[i]-mean, 2);
  std = sqrt(variance/newSize);

  // generate a probability histogram vector with width +/- 4*std
  vector<int> prob_histogram;
  // need to round the standard deviation and mean to the nearest +-4J
  double Emin = (mean - fmod(mean,4))-3*(std - fmod(std,4) + 4);
  // fmod replaces % for non-int values. Added four since the change is zero for std<4.
  double Emax = (mean - fmod(mean,4))+3*(std - fmod(std,4) + 4);
  // found no documentation on how to simplify this rounding to nearest 4th.
  // no solution to doubles rounding to nearest multiple of integer number

  prob_histogram.push_back(Emin);
  prob_histogram.push_back(Emax); // two first values tell about the domain.
  prob_histogram.push_back(mean);
  prob_histogram.push_back(std); // two next values tell about the distribution

  // loop the energy domain and count the number of times the energy E appears
  for (double Eval=Emin; Eval<=Emax; Eval+=4)
  {
    int counter=0;
    for (int i=startVal; i<total; i++) if (energy_vec[i]==Eval) counter++;
    prob_histogram.push_back(counter);
  }
  return prob_histogram;
}

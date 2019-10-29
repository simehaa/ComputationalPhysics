#include "lithosphere.h"

void JacobiMethodLithosphere(int situation){
  /*
   * Solver throughout 1 Ga of the lithosphere in three zones:
   * Zone 1: 0-20km,
   * Zone 2: 20-40km,
   * Zone 3: 40-120km.
   * Situation=1: no heat production
   * Situation=2: natural heat production
   * Situation=3: enriched mantle (zone 3) without decay
   * Situation=4: enriched mantle (zone 3) witho decay
   */

  string fn_base;
  string fn_end = ".bin";
  double mantleQ = 0; // represents the constant Q term in the mantle (zone 3)
  if (situation==1) {
    cout << "Solving lithosphere with no heat production\n";
    fn_base = "data/noheat/u";
  } else if (situation==2) {
    cout << "Solving lithosphere with natural heat production\n";
    fn_base = "data/naturalheat/u";
    mantleQ = 0.05;
  } else if (situation==3) {
    cout << "Solving lithosphere with enrichment without decay\n";
    fn_base = "data/nodecay/u";
    mantleQ = 0.55;
  } else if (situation==4) {
    cout << "Solving lithosphere with enrichment with decay\n";
    fn_base = "data/decay/u";
    mantleQ = 0.05;
  } else {
    cout << "Error, situation must be 1,2,3,4.\n";
    exit(1);
  }

  // initialization
  int nx = 125; // 0-150 km wide
  int ny = 100; // 0-120 km deep
  int nt = 100; // number of time steps
  double dx = 0.01; // step length [1.2 km]: 0.01*[1.2 km]*100 points = 120 km
  double dt = 0.01; // 1 = [10^9 yr]: 100*dt = 1 Ga
  mat u = zeros<mat>(nx+1,ny+1); // Scaled from 0-1 -> 8-1300 deg Celsius
  if (situation==1) {
    boundaryNoHeat(u,nx,ny);
  } else {
    boundaryNaturalHeat(u,nx,ny);
  }
  mat u_old = u;

  // Constants and variables
  int maxiter = 20000; // Max no. of iterations each time step in Jacobi Method
  double delta, diff, Q, time; // temporary constant that will be used
  double scale = (nx+1)*(ny+1); // No. of points in matrix
  double tol = 1e-10; // tolerance for convergence in Jacobi method
  double alpha = dt/(dx*dx);
  double rho = 3.51*1e3; // density [kg/m^3]
  double cp = 1000; // heat capacity [J/kg/deg Celsius]
  double k = 2.5; // thermal conductivity [W/m/deg Celsius]
  double uscale = 1292; // temperature scale
  double tscale = 3.1557e16; // time scale
  double xscale = 120000; // distance scale
  double beta = tscale*k/(rho*cp*xscale*xscale);
  double Qscale = 1e-6*tscale*dt/(rho*cp*uscale);
  double factor1 = beta*alpha;
  double factor2 = 1/(1 + 4*factor1);

  // Construction of Qvec: stores constant heat production in zones:
  vec Qvec = zeros<vec>(ny+1);
  if (situation != 1) {
    double Q1 = 1.4*Qscale;
    double Q2 = 0.35*Qscale;
    double Q3 = mantleQ*Qscale;
    for (int i=0; i<17; i++) Qvec(i) = Q1; // 0-20 km depth
    for (int i=17; i<34; i++) Qvec(i) = Q2; // 20-40 km depth
    for (int i=34; i<=ny; i++) Qvec(i) = Q3; // 40+ km depth
  }

  // EConstruction of Qtime: stores time dependet heat production (due to decay):
  vec Qtime = zeros<vec>(nt+1);
  if (situation==4) { // Qtime is simply 0 for cases 1-3
    double U_enrich, Th_enrich, K_enrich;
    for (int t=0; t<=nt; t++) {
      time = (double) t/nt;
      U_enrich = exp (-0.155*time);
      Th_enrich = exp (-0.0495*time);
      K_enrich = exp (-0.555*time);
      Qtime(t) = Qscale*(0.2*U_enrich + 0.2*Th_enrich + 0.1*K_enrich);
    }
  }

  // time loop
  for (int t=1; t<=nt; t++)
  {
    int iter = 0;
    double diff=1;
    mat u_guess = u; // Guess matrix, can be anything
    // Jacobi method loop
    while (iter < maxiter && diff > tol)
    {
      diff = 0;
      for (int j=1; j<ny; j++) {
          // Situation separation
          Q=Qvec(j); // constant Q for depth at index j
          if (j >= 34) Q += Qtime(t); // time dependent Q, only at depth j>=34
        for (int i=1; i<nx; i++) {
          // Jacobi method algorithm, u should converge towards solution
          delta = (u_guess(i,j+1)+u_guess(i,j-1)+u_guess(i+1,j)+u_guess(i-1,j));
          u(i,j) = factor2*(factor1*delta + Q + u_old(i,j));
          diff += fabs(u(i,j) - u_guess(i,j));
        }
      }
      u_guess = u;
      diff /= scale;
      iter++;
    } // end Jacobi method loop
    u_old = u;
    cout << "timestep: " <<  t << " iterations: " << iter << endl;
    u.save(fn_base+to_string(t)+fn_end, raw_binary); // Save matrix to proper folder
  } // end time loop
}

void boundaryNoHeat(mat& u, int nx, int ny) {
  // Creates boundary conditions for the case of Q = 0 everywhere
  for (int i=0; i<=nx; i++) {
    u(i,0) = 0;
    u(i,ny) = 1;
  }
  // linear temperature from 0-1 (scaled)
  for (int i=0; i<=nx; i+=nx) {
    for (double j=0; j<=ny; j++) {
      u(i,j) = j/ny;
      u(i,j) = j/ny;
    }
  }
}

void boundaryNaturalHeat(mat& u, int nx, int ny) {
  // Analytic solution to steady state Temp(depth) with natural heat production
  // temperature is scaled [8-1300] -> [0-1292] -> [0,1]
  double y, temp;
  for (double j=0; j<17; j++) {
    y = j*1.2;
    temp = (-0.28*y*y + 23.66*y)/1292.0;
    for (int i=0; i<=nx; i++) u(i,j) = temp;
  }
  for (double j=17; j<34; j++) {
    y = j*1.2;
    temp = (-0.07*y*y + 15.26*y + 86)/1292.0;
    for (int i=0; i<=nx; i++) u(i,j) = temp;
  }
  for (double j=34; j<=ny; j++) {
    y = j*1.2;
    temp = (-0.01*y*y + 10.46*y + 180)/1292.0;
    for (int i=0; i<=nx; i++) u(i,j) = temp;
  }
}

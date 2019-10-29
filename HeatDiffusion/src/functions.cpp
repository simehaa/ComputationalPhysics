#include "functions.h"

void analytic1D(int nx, double dx, string filename) {
  /*
   * Analytic solution for the one dimensional diffusion equation,
   * with boundary conditions:
   * u(L,t) = 1 and u(0,t) = 0, for all t.
   * u(x,0) = 0 for x < L.
   */
  int infty = 1000; // what is defined as numerical "infinity".
  float L=1; // scale the rod such that x goes from 0 to L=1.
  mat u = zeros<mat>(nx+1, 2);
  double an, lambdan, lambdanL, factor;
  int c = 0;
  // time loop
  for (double t=0.1; t<=0.2; t+=0.1){
    // position x loop
    for (int x = 0; x<=nx; x++){
      // calculate the transient solution.
      factor = 0;
      for (int n = 1; n<=infty; n++){
        lambdanL = M_PI*n;
        lambdan = lambdanL/L;
        an = 2*(lambdanL*cos(lambdan*L) - sin(lambdanL))/(lambdan*lambdanL);
        factor += an * exp(-lambdan*lambdan*t) * sin(lambdan*x*dx);
      }
      u(x,c) = x*dx/L + factor;
    }
    c++; // Easter egg
  }
  u.save(filename, raw_binary);
}

void analytic2D(int nx, double dx, double t, double dt, string filename){
  /*
   * Analytic solution for the two dimensional diffusion equation at a
   * given time. Boundary conditions are all equal to zero. A "sine-paraboloid"
   * which is centered in the middle of the medium is the initial condition.
   * Boundary conditions are zeros.
   */
  mat u = zeros<mat>(nx+1, nx+1); // Assume a square grid
  for (int x=0; x<=nx; x++){
    for (int y=0; y<=nx; y++){
      u(x,y) = exp(-2*M_PI*M_PI * t*dt) * sin(M_PI*x*dx) * sin(M_PI*y*dx);
    }
  }
  u.save(filename, raw_binary);
  cout << "File " << filename << " written." << endl;
}

void explicitForwardEuler(int nx, int nt, double dx, double dt, string filename) {
  /*
   * Solves position in 1D for all times
   * using the explicit forward Euler scheme: straight forward method which
   * is simply solved in a position loop in a time loop/
   */

  // initialization
  int i,t,tp;
  double alpha = dt/dx/dx;
  if (alpha > 0.5) {
    cout << "Warning: alpha too large for explicit scheme \n";
  }
  double beta = (1 - 2*alpha);
  mat u = zeros<mat>(nx+1,nt+1);

  // boundary conditions, [u(x,0)=0 and u(0,t)=0]
  for (t=0; t<=nt; t++) u(nx,t) = 1;

  // time loop
  for (t=1; t<=nt; t++) {
    tp = t-1; // previous time step
    for (i=1; i<nx; i++) { // inner time points
      u(i,t) = alpha*(u(i-1,tp) + u(i+1,tp)) + beta*u(i,tp);
    }
  }
  u.save(filename, raw_binary);
}

void implicitBackwardEuler(int nx, int nt, double dx, double dt, string filename) {
  /*
   * Solves position in 1D for all times using
   * the implicit backward Euler scheme. Calls tridiagonalSolver with
   * beta as diagonal elements, and -alpha as offdiagonal elements.
   */

  // initialization
  int i,t,tn;
  double alpha = dt/dx/dx;
  double beta = (1 + 2*alpha); // diagonal elements
  alpha *= -1; // offdiagonal elements
  mat u = zeros<mat>(nx+1,nt+1);
  vec u_new = zeros<vec>(nx+1);

  // boundary conditions, [u(x,0)=0 and u(0,t)=0]
  for (t=0; t<=nt; t++) u(nx,t) = 1;

  // time loop
  for (t=1; t<=nt; t++) {
    u_new = u.col(t);
    tridiagonalSolver(u_new, u.col(t-1), beta, alpha, nx);
    u.col(t) = u_new;
  }
  u.save(filename, raw_binary);
}

void CrankNicolsonScheme(int nx, int nt, double dx, double dt, string filename) {
  /*
   * Uses a combination of the explicit and implicit scheme: Calls function
   * tridiagonalSolver with a different tridiagonal matrix than the implicit
   * scheme, this matrix is derived by hand.
   */

  // initialization
  int i,t,tn,tp;
  double alpha = dt/dx/dx;
  double beta = (2 - 2*alpha);
  double diag = (2 + 2*alpha);
  double offdiag = -alpha;
  mat u = zeros<mat>(nx+1,nt+1);
  vec b = zeros<vec>(nx+1);
  vec u_new = zeros<vec>(nx+1);

  // boundary conditions
  for (t=0; t<=nt; t++) u(nx,t) = 1;

  // time loop
  for (t=1; t<=nt; t++) {
    tp = t-1; // previous time step

    // Explicit part
    for (i=1; i<nx; i++) { // inner time points
      b(i) = alpha*(u(i-1,tp) + u(i+1,tp)) + beta*u(i,tp);
    }
    // Implicit part
    u_new = u.col(t);
    tridiagonalSolver(u_new, b, diag, offdiag, nx);
    u.col(t) = u_new;
  }
  u.save(filename, raw_binary);
}

void tridiagonalSolver(vec& u, vec b, double diag, double offdiag, int n) {
  /*
   * Thomas algorithm:
   * Solves matrix vector equation Au = b,
   * for A being a tridiagonal matrix with constant
   * elements diag on main diagonal and offdiag on the off diagonals.
   */
  vec beta = zeros<vec>(n+1); beta[1] = diag;
  vec u_old = zeros<vec>(n+1); u_old(1) = b(1);
  double btemp;

  // forward substitution
  for(int i=2; i<=n; i++){
    btemp = offdiag/beta[i-1];
    beta(i) = diag - offdiag*btemp;
    u_old(i) = b(i) - u_old(i-1)*btemp;
  }

  // Special case, boundary conditions
  u(0) = 0;
  u(n) = 1;

  // backward substitution
  for(int i=n-1; i>0; i--){
    u(i) = (u_old(i) - offdiag*u(i+1))/beta(i);
  }
}

void JacobiMethod(){
  /*
   * Iterative solver for a of a diagonally dominant system og linear equations.
   */
  // initialization
  int nx = 1000;
  int nt = 20;
  double dx = 0.001;
  double dt = 0.0001;

  // Initial conditions
  mat u = zeros<mat>(nx+1,nx+1);
  for (double i=1; i<nx; i++) {
    for (double j=1; j<nx; j++) {
      u(i,j) = sin(i*M_PI/(double)nx)*sin(j*M_PI/(double)nx);
    }
  }

  // Constants, variables
  mat u_old = u;
  int maxiter = 1000;
  double delta, tol=1e-8;
  double alpha = dt/(dx*dx);
  double factor = 1.0/(1.0 + 4*alpha);
  double factor_a = alpha*factor;
  double scale = (nx+1)*(nx+1);
  string fn_base = "data/twodim_dx0001_dt00001/u";
  string filename;
  string fn_end = ".bin";

  // time loop
  for (int t=1; t<=nt; t++)
  {
    int iter = 0;
    double diff=1;
    mat u_guess = ones<mat>(nx+1,nx+1);
    while (iter < maxiter && diff > tol)
    {
      diff = 0;
      // Loop over all inner elements, will converge towards solution
      for (int j=1; j<nx; j++) {
        for (int i=1; i<nx; i++) {
          // u_guess is the previous u, which also work for a random guess
          delta = (u_guess(i,j+1)+u_guess(i,j-1)+u_guess(i+1,j)+u_guess(i-1,j));
          u(i,j) = factor_a*delta + factor*u_old(i,j);
          diff += fabs(u(i,j) - u_guess(i,j));
        }
      } // end of double for loop
      u_guess = u;
      diff /= scale;
      iter++;
    } // end iteration loop
    filename = fn_base + to_string(t) + fn_end;
    u_old = u;
    u.save(filename, raw_binary);
    cout << t << endl;
  } // end time loop
  ofstream ofile;
  ofile.open("data/twodimensions.txt");
  ofile << nx << endl;
  ofile << nt << endl;
  ofile.close();
}

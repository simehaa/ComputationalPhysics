#include "functions.h"
#include <cmath>
#include <math.h>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <time.h>

using namespace std;
using namespace arma;


void analytic_eigenvalues(mat A){
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, A); // eigenvalues for a symmetric matrix
  eigval.print();
}

double max_value_indexes(mat A, int N, int& k, int& l){
  double maxval=0, val;
  for (int i=0; i<(N-1); i++){ // nested loop over upper triangle of A (symmetry)
    for (int j=i+1; j<(N-1); j++){
        val = A(i, j);
        if (fabs(val)>maxval){ // test for largest absolute value
          k=i; l=j; // A(k,l) has the largest absolute value
          maxval = fabs(val);
        }
      }
    }
  return maxval;
}

mat generate_A_matrix(int N, vec a, vec d){
  mat A = zeros<mat>(N-1,N-1);
  A(0,0) = d(0); A(0,1) = a(0); // first row
  for (int i=1; i<N-2; i++) { // All inner rows in tridiag matrix filled row by row
    A(i, i) = d(i);
    A(i,i+1) = a(i);
    A(i,i-1) = a(i);
  }
  A(N-2,N-2) = d(N-2); A(N-2,N-3) = a(N-3); // last row
  return A;
}

void Jacobi_Rotation_algorithm(mat& A, mat& R, int N, int k, int l){
  // 35 FLOPS
  // Obtaining values tau, t (tan), c (cos), s (sin)
  double tau, t, s, c, a_ik, a_il, a_kk, a_ll, r_ik, r_il;
  tau = (A(l,l) - A(k,k))/(2*A(k,l));
  if (tau >= 0) {
    t = 1.0/(tau + sqrt(1 + tau*tau));
  } else {
    t = -1.0/(-tau + sqrt(1 + tau*tau));
  }
  c = 1/(sqrt(1 + t*t));
  s = c*t;
  // reduces three flops in total:
  double cc = c*c;
  double ss = s*s;
  double cs = c*s;

  a_kk = A(k,k)*cc - 2*A(k,l)*cs + A(l,l)*ss;
  A(l,l) = A(l,l)*cc + 2*A(k,l)*cs + A(k,k)*ss;
  A(k,l) = 0; // required
  A(l,k) = 0;
  A(k,k) = a_kk;

  for (int i=0; i<N-1; i++) {
    if (i != k && i != l) {
      a_ik = A(i,k)*c - A(i,l)*s;
      a_il = A(i,l)*c + A(i,k)*s;
      A(i,k) = a_ik;
      A(i,l) = a_il;
      A(k,i) = a_ik;
      A(l,i) = a_il;
    }
    r_ik = R(i,k);
    r_il = R(i,l);
    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
  }
}

void find_lowest_eigval_eigvec_pair(double& eigval, vec& eigvec, mat A, mat A_original, mat R, int N) {
  int minIndex = 0;
  for(int i=1; i<N-1; i++){
    // find the index of the smallest element
    if ( A(i,i) < A(minIndex,minIndex) ) {
      minIndex = i;
      eigval = A(i,i);
    }
  }
  // test that eigenpair is correct:
  eigvec = R.col(minIndex);
  vec x = A_original*eigvec;
  vec y = eigval*eigvec;
  double tolerance = 1e-10;
  for (int i=0; i<N-1; i++) { // x and y should analytically be the same
    if (fabs(x(i) - y(i)) > tolerance) {
      cout << "Error: eigenvalue/eigenvector pair is incorrect (A*x != lambda*x), with tolerance: " << tolerance << endl;
      cout << "index: " << i << " value: " << y(i) << " error: " << fabs(x(i) - y(i)) << endl;
      exit(1);
    }
  }
}

void write_file(int N, int j, double eigval, double wr, vec rho, vec eigvec) {
  string filename = "./data/project2_";
  filename.append(to_string(j) + ".txt");
  ofstream ofile;
  ofile.open(filename, std::ofstream::out | std::ofstream::trunc);

  ofile << setw(10) << "lambda: " << eigval << endl;
  ofile << setw(10) << "wr: " << wr << endl;
  ofile << setw(20) << "rho:" << setw(20) << "eigvec: " << endl;

  for (int i=0; i<N-1; i++) {
    ofile << setw(20) << setprecision(10) << rho(i);
    ofile << setw(20) << setprecision(10) << eigvec(i) << endl;
  }
  ofile.close();
}

void buckling_beam(int N) {
  double h = 1.0/N;
  vec a = ones<vec>(N-2); // upper and lower diagonals
  vec d = ones<vec>(N-1); // main diagonal
  a *= -1.0/(h*h);
  d *= 2.0/(h*h);
  mat A = generate_A_matrix(N, a, d);
  mat R = zeros<mat>(N-1,N-1);

  // Jacobi's algorithm
  double maxvalue = 10.;
  double epsilon = 1e-16;
  int iteration = 0;
  int maxIterations = 100000;
  int k, l;
  // loop until nondiagonal maxvalue is smaller than epsilon OR max iterations
  double start, finish;
  start = clock();
  while ( maxvalue > epsilon && iteration < maxIterations ) { // Main algorithm loop that performs rmatrix otations
    maxvalue = max_value_indexes(A, N, k, l);
    iteration++;
    Jacobi_Rotation_algorithm(A, R, N, k, l);
  }
  finish = clock();
  double timeElapsed = (finish-start)/CLOCKS_PER_SEC;

  // storing and sorting eigenvalues:
  vec eigvals = zeros<vec>(N-1);
  for (int i=0; i<N-1; i++) {
    eigvals(i) = A(i,i);
  }
  sort(eigvals.begin(), eigvals.end());

  // terminal print:
  cout << "  time: " << timeElapsed << " s" << endl;
  cout << "  iterations: " << iteration << endl;
  cout << "\n  numerical eigvals: analytic eigvals:" << endl;
  for (int j=0; j<N-1; j++) {
    cout << "  " << setw(17) << setprecision(16) << left << eigvals(j) << "  ";
    cout << d(0) + 2*a(0)*cos((j+1)*M_PI/N) << endl;
  }
}

void one_electron_system(int N, double rho_max) {
  double rho_0 = 0;
  double h = (rho_max - rho_0)/N;
  vec a = ones<vec>(N-2); // upper and lower diagonals
  vec d = ones<vec>(N-1); // main diagonal
  a *= -1.0/(h*h);
  d *= 2.0/(h*h);
  double rho_i;
  for (int i=0; i<N-1; i++) {
    rho_i = rho_0 + (i+1)*h;
    d(i) += rho_i*rho_i;
  }

  double maxvalue = 10.;
  double epsilon = 1e-8;
  int iteration = 0;
  int maxIterations = 2000000;
  int k, l;
  mat R = eye<mat>(N-1,N-1); // matrix where columns will store eigenvectors
  mat A = generate_A_matrix(N, a, d); // this will be changed
  // loop until nondiagonal maxvalue is smaller than epsilon OR max iterations
  double start, finish;
  start = clock();
  while ( maxvalue > epsilon && iteration < maxIterations ) { // Main algorithm loop that performs rmatrix otations
    maxvalue = max_value_indexes(A, N, k, l);
    iteration++;
    Jacobi_Rotation_algorithm(A, R, N, k, l);
  }
  finish = clock();
  double timeElapsed = (finish-start)/CLOCKS_PER_SEC;

  vec eigvals = zeros<vec>(N-1);
  for (int i=0; i<N-1; i++) {
    eigvals(i) = A(i,i);
  }
  sort(eigvals.begin(), eigvals.end());
  // terminal print:
  cout << "  time: " << timeElapsed << " s" << endl;
  cout << "  iterations: " << iteration << endl;
  cout << "\n  first 5 eigenvalues:" << endl;
  for (int j=0; j<5; j++) {
    cout << "  " << eigvals(j) << endl;
  }
}

void two_electron_system(int N, double rho_max) {
  vec wrvec = zeros<vec>(4);
  wrvec(0) = 0.01;
  wrvec(1) = 0.5;
  wrvec(2) = 1.0;
  wrvec(3) = 5.0;

  double rho_0 = 0;
  double h = (rho_max - rho_0)/N;
  vec a = ones<vec>(N-2); // upper and lower diagonals
  a *= -1.0/(h*h);

  // loop over different omega_r = wr
  for (int j=0; j<4; j++) {
    double wr = wrvec(j);
    vec d = ones<vec>(N-1); // main diagonal
    d *= 2.0/(h*h);

    // adding of the extra term V (potential) to the main diagonal d
    vec rho = zeros<vec>(N-1);
    for (int i=0; i<N-1; i++) { // loop over diagonal elements:
      rho(i) = rho_0 + (i+1)*h;
      d(i) += wr*wr*rho(i)*rho(i) + 1.0/rho(i); // potential term
    }

    // Jacobi's method:
    mat R = eye<mat>(N-1,N-1); // matrix where columns will store eigenvectors
    mat A = generate_A_matrix(N, a, d); // this will be changed
    mat A_original = generate_A_matrix(N, a, d); // this will remain unchanged
    double maxvalue = 10.;
    double epsilon = 1e-10;
    int iteration = 0;
    int maxIterations = 100000;
    int k, l;
    // loop until nondiagonal maxvalue is smaller than epsilon OR max iterations
    double start, finish;
    start = clock();
    while ( maxvalue > epsilon && iteration < maxIterations ) { // Main algorithm loop that performs rmatrix otations
      maxvalue = max_value_indexes(A, N, k, l);
      iteration++;
      Jacobi_Rotation_algorithm(A, R, N, k, l);
    }
    finish = clock();
    double timeElapsed = (finish-start)/CLOCKS_PER_SEC;
    cout << "\n  Jacobi's method done, number of iterations: " << iteration << endl;
    cout << "  time: " << timeElapsed << " s" << endl;

    // finding the ground state eigenpair:
    double eigval;
    vec eigvec = zeros<vec>(N-1);
    find_lowest_eigval_eigvec_pair(eigval, eigvec, A, A_original, R, N);
    cout << "  lowest eigenpair found, with eigenvalue: " << eigval << endl;

    // file writing
    write_file(N, j, eigval, wr, rho, eigvec);
    cout << "  file written." << endl;
    cout << "  wr = " << wr << " is done." << endl;
  }
}

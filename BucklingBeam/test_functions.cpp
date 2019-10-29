#include <armadillo>
#include <iostream>
#include <vector>
#include <math.h>
#include "test_functions.h"
#include "functions.h"

using namespace std;
using namespace arma;

void test_max_value_indices() {
  /* Generating a 9x9 matrix, and
  one element is set to -5.0, with
  known indiced */
  int N = 10;
  vec a = ones<vec>(N-2) * -1.0;
  vec d = ones<vec>(N-1) * 2.0;
  mat A = generate_A_matrix(N, a, d);
  A(2,4) = -5.0; // indices k = 2, l = 4
  int k,l;
  double maxval;
  maxval = max_value_indexes(A,N,k,l);
  if (maxval != 5.0 || k != 2 || l != 4) {
    cout << "test_max_value_indexes failed!" << endl;
    exit(1);
  }
}

void test_eigenvalues() {
  /* Generating a 3x3 matrix,
  then performing jacobi's algo,
  then sorting the eigenvalues from smallest to largest.
  finally testing the eigenvalues vs. analytic eigenvalues */
  int N=4;
  int explode=50000;
  int iteration=0;
  int k, l;
  vec a = ones<vec>(N-2) * -1.0;
  vec d = ones<vec>(N-1) * 2.0;
  mat A = generate_A_matrix(N, a, d); // has analytic eigenvalues
  mat R = eye<mat>(N-1,N-1);
  double maxvalue=10.0;
  double epsilon=1e-12;
  // Jacobi algorithm
  while(maxvalue>epsilon && iteration<=explode){
    maxvalue = max_value_indexes(A, N, k, l);
    iteration++;
    Jacobi_Rotation_algorithm(A, R, N, k, l);
  }
  // Test of analytical vs. numerical eigenvalues
  vector<double> numerical_eigenvalues;
  for (int i=0; i<3; i++) {
    numerical_eigenvalues.push_back(A(i,i));
  }
  // three eigenvealues sorted from smallest to greatest
  sort(numerical_eigenvalues.begin(), numerical_eigenvalues.end());
  // test of all three eigenvalues
  for (int i=1; i<N; i++) {
    double ana_eigval = 2.0 + 2.0*-1.0*cos(i*M_PI/N); // divided by N for an (N-1) x (N-1) matrix
    double num_eigval = numerical_eigenvalues[i-1];
    if (fabs(num_eigval - ana_eigval) > 1e-8) { // test with a tolerance 1e-8
      cout << "test_eigenvalues failed: ";
      cout << num_eigval << " != ";
      cout << ana_eigval << endl;
      exit(1);
    }
  }
}

void test_orthogonality() {
  int N=6;
  int iteration=0;
  int k, l;
  vec a = ones<vec>(N-2) * -1.0;
  vec d = ones<vec>(N-1) * 2.0;
  mat A = generate_A_matrix(N, a, d); // has analytic eigenvalues
  double maxvalue=10.0;
  double epsilon=1e-12;
  mat R = eye<mat>(N-1,N-1);
  // Jacobi algorithm for 70 iterations
  for (int j=0; j<70; j++){
    maxvalue = max_value_indexes(A, N, k, l);
    iteration++;
    Jacobi_Rotation_algorithm(A, R, N, k, l);
  }
  // test of dot product of column 2 and 5
  double dotProduct;
  for (int i=0; i<N-1; i++) {
    dotProduct += A(i,1)*A(i,4);
  }
  if (fabs(dotProduct) > 1e-8) {
    cout << "test_orthogonality failed: 2nd and 5th";
    cout << " column were not orthogonal after 70 iterations" << endl;
    exit(1);
  }
}

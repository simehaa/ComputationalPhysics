#include <iostream>
#include "test_functions.h"
#include "functions.h"

using namespace std;

int main(int argc, char* argv[])
{
  // test functions:
  test_max_value_indices();
  test_eigenvalues();
  test_orthogonality();
  cout << "Test functions passed." << endl;

  // the three different systems:
  cout << "\nBuckling beam system:" << endl;
  //buckling_beam(10);

  cout << "\nOne electron in a harmonic oscillator:" << endl;
  //one_electron_system(200, 10.0);

  cout << "\nTwo electrons in a harmonic oscillator:" << endl;
  two_electron_system(200, 10.0);

  cout << "\nDone." << endl;
  return 0;
}

// Terminal run:

/*
simen@simen-ubuntu:~/steinngithub/FYS3150/Project2$ make -j
g++ -c functions.cpp
g++ -o test.x main.o functions.o test_functions.o -llapack -larmadillo -lblas
simen@simen-ubuntu:~/steinngithub/FYS3150/Project2$ ./test.x
Test functions passed.

Buckling beam system:
  time: 0.000242 s
  iterations: 143

  num. eigval: ana. eigval:
  9.7887       9.7887
  38.1966      38.1966
  82.4429      82.4429
  138.197      138.197
  200          200
  261.803      261.803
  317.557      317.557
  361.803      361.803
  390.211      390.211

One electron in a harmonic oscillator:
  time: 13.8185 s
  iterations: 63410

  first 5 eigenvalues:
  2.99922
  6.99609
  10.9905
  14.9823
  18.9717

Two electrons in a harmonic oscillator:

  Jacobi's method done, number of iterations: 70691
  time: 15.4065 s
  lowest eigenpair found, with eigenvalue: 0.840745
  file written.
  wr = 0.01 is done.

  Jacobi's method done, number of iterations: 70949
  time: 15.6242 s
  lowest eigenpair found, with eigenvalue: 2.23097
  file written.
  wr = 0.5 is done.

  Jacobi's method done, number of iterations: 70536
  time: 15.6938 s
  lowest eigenpair found, with eigenvalue: 4.05767
  file written.
  wr = 1 is done.

  Jacobi's method done, number of iterations: 69056
  time: 15.165 s
  lowest eigenpair found, with eigenvalue: 17.4436
  file written.
  wr = 5 is done.

Done.
*/

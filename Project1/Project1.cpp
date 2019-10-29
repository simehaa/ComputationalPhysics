#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <armadillo>
using namespace std;
using namespace arma;
ofstream ofile;
// test
void ThomasAlgorithmGeneralized(string fname, int maxexp);
void LUDecompFunc(string fname, int maxexp);
void ThomasAlgorithmSpecialized(string fname, int maxexp);

int main(int argc, char* argv[])
{
    string fname;
    int maxexp;
    fname=argv[1];
    maxexp = atoi(argv[2]);

    //ThomasAlgorithmGeneralized(fname, maxexp);
    //LUDecompFunc(fname, maxexp);
    ThomasAlgorithmSpecialized(fname, maxexp);
}


void ThomasAlgorithmGeneralized(string fname, int maxexp){
  clock_t start, finish;
  string filename;
  for(int i=1; i<=maxexp; i++){
    filename=fname;
    int n;
    n=pow(10.0, i);
    filename.append(to_string(i));
    double *d_vec = new double [n+2]; double *dt_vec = new double [n+2]; double *f_vec = new double [n+2];
    double *ft_vec = new double [n+2]; double *a_vec = new double [n+1]; double *c_vec = new double [n+1]; double *u_vec = new double [n+2];

    double x0=0, h, x;
    h=1/(double(n+1));

    //Fill f array analytically:
    for (int i=0; i<=n+1; i++){
        x=x0+h*double(i);
        f_vec[i]=h*h*100*exp(-10*x);
    }
    //Special case diagonals
    for(int i=1; i<=n; i++){
        d_vec[i]=2;
        a_vec[i]=-1;
        c_vec[i]=-1;
    }

    dt_vec[1]=d_vec[1];
    ft_vec[1]=f_vec[1];
    u_vec[0]=0;
    u_vec[n+1]=0;

    start = clock();

    //First for-loop to create b vector.
    for(int i=2; i<=n; i++){
        dt_vec[i]=d_vec[i] - c_vec[i-1]*a_vec[i-1]/dt_vec[i-1];
        ft_vec[i]=f_vec[i] - ft_vec[i-1]*a_vec[i-1]/dt_vec[i-1];
    }

    u_vec[n]=ft_vec[n]/dt_vec[n];

    for(int i=n; i>0; i--){
        u_vec[i]=(ft_vec[i]-c_vec[i]*u_vec[i+1])/dt_vec[i];
    }

    finish = clock();
    double clocks = (finish-start);
    double timeElapsed = clocks/CLOCKS_PER_SEC;
    cout << "N=10^" << i << ": " << timeElapsed << "s" << endl;


    //write file with current i and n:
    filename.append(".txt");
    ofile.open(filename);
    //ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(double i=1; i<=n; i++){
        double h;
        h=1/(double(n)+1);
        ofile << setprecision(10) << setw(20) << h*i;
        ofile << setprecision(10) << setw(20) << 1-(1-exp(-10))*h*i - exp(-10*h*i);
        ofile << setprecision(10) << setw(20) << u_vec[int(i)];
        ofile << setprecision(10) << setw(20) << log10(fabs((- (1-(1-exp(-10))*h*i - exp(-10*h*i)) + u_vec[int(i)])/(1-(1-exp(-10))*h*i - exp(-10*h*i)))) << endl;
    }
    ofile.close();

    delete [] d_vec; delete [] dt_vec; delete [] f_vec; delete [] ft_vec;
    delete [] a_vec; delete [] c_vec; delete [] u_vec;
  }
  cout << "Calculation complete." << endl;
}

void LUDecompFunc(string fname, int maxexp){
  clock_t start, finish;
  string filename;
  for (int i=1; i<=maxexp; i++){
    filename = fname;
    filename.append(to_string(i));
    int n=pow(10.0,i);
    double h;
    h=1./(n+1);
    mat A = zeros<mat>(n,n);

    for(int i=0; i<n; i++){
      A(i,i)=2;
    }
    for(int i=0; i<(n-1); i++){
      A(i,i+1)=-1;
      A(i+1,i)=-1;
    }

    vec b = zeros<vec>(n);
    for(int i=0; i<n; i++){
      b(i)=100*exp(-10*i*h)*h*h;
    }

    start = clock(); //Algorithm start
    mat L, U;
    vec x, z, u;
    lu(L,U,A);

    z = solve(L, b);
    u = solve(U, z);
    finish = clock(); //Algorithm completed

    double clocks = (finish-start);
    double timeElapsed = clocks/CLOCKS_PER_SEC;
    cout << "N=10^" << i << ": " << timeElapsed << "s" << endl;


    filename.append(".txt");
    ofile.open(filename);
    //ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(double i=1; i<n; i++){
        double h;
        h=1/(double(n)+1);
        ofile << setprecision(10) << setw(20) << h*i;
        ofile << setprecision(10) << setw(20) << 1-(1-exp(-10))*h*i - exp(-10*h*i);
        ofile << setprecision(10) << setw(20) << u(int(i));
        ofile << setprecision(10) << setw(20) << log10(fabs((- (1-(1-exp(-10))*h*i - exp(-10*h*i)) + u(int(i)))/(1-(1-exp(-10))*h*i - exp(-10*h*i)))) << endl;
    }
    ofile.close();

  }
}

void ThomasAlgorithmSpecialized(string fname, int maxexp){
  clock_t start, finish;
  string filename;
  for(int i=1; i<=maxexp; i++){
    filename=fname;
    int n;
    n=pow(10.0, i);
    filename.append(to_string(i));
    double *dt_vec = new double [n+2]; double *f_vec = new double [n+2];
    double *ft_vec = new double [n+2]; double *u_vec = new double [n+2];

    double x0=0, h, x;
    h=1/(double(n+1));

    for (int i=0; i<=n+1; i++){
        x=x0+h*double(i);
        f_vec[i]=h*h*100*exp(-10*x);
    }
    dt_vec[1]=2;
    ft_vec[1]=f_vec[1];
    u_vec[0]=0;
    u_vec[n+1]=0;

    start = clock();
    for(int i=2; i<=n; i++){
        dt_vec[i]=2 - 1/dt_vec[i-1];
        ft_vec[i]=f_vec[i] + ft_vec[i-1]/dt_vec[i-1];
    }
    u_vec[n]=ft_vec[n]/dt_vec[n];

    for(int i=n; i>0; i--){
        u_vec[i]=(ft_vec[i]+u_vec[i+1])/dt_vec[i];
    }

    finish = clock();
    double clocks = (finish-start);
    double timeElapsed = clocks/CLOCKS_PER_SEC;
    cout << "N=10^" << i << ": " << timeElapsed << "s" << endl;


    //write file with current i and n:
    filename.append(".txt");
    ofile.open(filename);
    //ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(double i=1; i<=n; i++){
        double h;
        h=1/(double(n)+1);
        ofile << setprecision(10) << setw(20) << h*i;
        ofile << setprecision(10) << setw(20) << 1-(1-exp(-10))*h*i - exp(-10*h*i);
        ofile << setprecision(10) << setw(20) << u_vec[int(i)];
        ofile << setprecision(10) << setw(20) << log10(fabs((- (1-(1-exp(-10))*h*i - exp(-10*h*i)) + u_vec[int(i)])/(1-(1-exp(-10))*h*i - exp(-10*h*i)))) << endl;
    }
    ofile.close();

    delete [] dt_vec; delete [] f_vec; delete [] ft_vec;
    delete [] u_vec;
  }
  cout << "Calculation complete." << endl;
}


/*
Generalized Thomas alogrithm:

steinn@steinn-VirtualBox:~/Desktop/FYS3150/Project1$ ./a.out results 7
N=10^1: 4e-06s
N=10^2: 4e-06s
N=10^3: 3.1e-05s
N=10^4: 0.000327s
N=10^5: 0.003473s
N=10^6: 0.038868s
N=10^7: 0.354114s
Calculation complete.


LU - decomposition:

steinn@steinn-VirtualBox:~/Desktop/FYS3150/Project1$ ./a.out results 4
N=10^1: 0.000117s
N=10^2: 0.000811s
N=10^3: 0.290671s
N=10^4: 307.447s


Specialized Thomas algorithm:

steinn@steinn-VirtualBox:~/Desktop/FYS3150/Project1$ ./a.out results 7
N=10^1: 3e-06s
N=10^2: 4e-06s
N=10^3: 2.7e-05s
N=10^4: 0.000316s
N=10^5: 0.003318s
N=10^6: 0.035397s
N=10^7: 0.334385s
Calculation complete.

*/

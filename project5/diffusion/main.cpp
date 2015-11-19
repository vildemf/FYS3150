#include <iostream>
#include <armadillo>
#include <time.h>
#include "lib.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace arma;

void forward_euler(int, int, double, vec&, vec&, vec&);
double initial(double);

int main()
{

    ofstream outfile;
    outfile.open("/home/vilde/Documents/FYS3150/project5/test.txt");
    double T = 0.5; // t(Nt)
    double d = 1;  // x(Nx)
    int Nt = 50000;   // number of t-values
    int Nx = 200;   // number of x-values
    vec t = linspace<vec>(0, T, Nt+1);
    vec x = linspace<vec>(0, d, Nx+1);
    double dt = t(1) - t(0);
    double dx = x(1) - x(0);

    double D = 1;
    double alpha = D*dt/(dx*dx);

    vec v_prev = zeros<vec>(Nx+1); // store previous for computing for next t
    vec v = zeros<vec>(Nx+1);

    forward_euler(Nt, Nx, alpha, v_prev, v, x);

    /*
     *
     */
    outfile << Nt << " " << Nx << " " << dt << " " << dx << " " << T << endl;
    for (int i=0; i<=Nx; i++) {
        outfile << setprecision(8) << v(i) << endl;
    }




    cout << "Hello World!" << endl;
    return 0;
}

/*
 * Stability condition forward_euler: dt/(dx**2) < 0.5
 * gives the condition Nt > 2T(Nx**2)
 *
 */
void forward_euler(int Nt, int Nx, double alpha, vec& v_prev, vec& v, vec& x) {
    // t=0
    // v(0), v(Nx) = 0 or initial() at t=0?
    for (int i=1; i<Nx; i++) {
        v_prev(i) = initial(x(i));
    }


    // Calculate v(x, t) for t=dt,...,T and x=dx,..,d-dx
    // Corresponds to j=1,...,Nt and i=1,...,Nx-1
    // Leave out i=0 and i=Nx to keep v here 0 (boundary condition)
    for (int j=1; j<=Nt; j++) {
        for (int i=1; i<Nx; i++) {
            v(i) = alpha*v_prev(i-1) + (1-2*alpha)*v_prev(i) + alpha*v_prev(i+1);
        }
        for (int i=0; i<=Nx; i++) {
            v_prev(i) = v(i);
        }
    }
}


double initial(double x) {
    return x-1;
}

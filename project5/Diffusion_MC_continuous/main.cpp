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

void mc_sampling(int, double, int, double, double, vec&, double, vec&);
// function for gaussian random numbers
double gaussian_deviate(long *);
void output(int, int, vec&, vec&);
void set_up(char**, double&, double&, int&);


/* Main:
 * This program solves the diffusion equation by Monte Carlo simulations and random walk.
 * The position values are continuous.
 * Command line input (in correct order):
 * outfilename, T, dt, mc_number
 */
int main(int argc, char *argv[])
{
    int mc_number = 0;
    double T, dt=0;
    set_up(argv, T, dt, mc_number);
    double move_prob=0.5;
    double D=1;
    double d=1;
    int Nx = 100;
    int Nt = round(T/dt) + 1;
    vec x = linspace<vec>(0, d, Nx+1);
    double dx = x(1) - x(0);

    // Array to represent the concentration
    vec pos_count = zeros<vec>(Nx+1);

    // Do the mc sampling
    mc_sampling(mc_number, T, Nx, dt, move_prob, pos_count, D, x);

    // Write to file
    ofstream ofile(argv[1]);
    for (int i=0; i<=Nx; i++) {
        ofile << pos_count(i)/(pos_count(0)) << endl; // Normalize
    }
    ofile.close();

    return 0;
}

void set_up(char** argv, double& T, double& dt, int& mc_number) {
    T = atof(argv[2]);
    dt = atof(argv[3]);
    mc_number = atof(argv[4]);

}

void mc_sampling(int mc_number, double T, int Nx, double dt, double move_prob, vec& pos_count, double D, vec& x) {
    // Declerations
    double l0 =0;
    double ksi = 0;
    double pos = 0;
    int pos_index = 0;
    long idum1 = -1;
    long idum2 = -1;
    for (int particles=1; particles<=mc_number; particles++) {
        // Particle's position when it enters the cleft
        pos_index = 0;
        pos = ran0(&idum1)*(x(1)-x(0));
        pos_count(pos_index) += 1;

        /*
         * Count the particle's position at every time step.
         * This corresponds to storing the movement of particles that had a time interval smaller than T
         * This corresponds to incorporating the number of particles entering the system at a time later than t=0,
         * having less than time T in the cleft
         */
        for (double time = dt; time <= T; time+=dt) {
            // new step
            ksi = gaussian_deviate(&idum2);
            l0 = ksi*sqrt(2*D*dt);
            pos += l0;

            // if it has reached the end, stop counting it
            if (pos >= x[Nx]) break;

            // if above interval, update index
            if (pos_index!=Nx) {
                if (pos > x[pos_index + 1]) {
                    pos_index += 1;
                }
            }

            // if it has reached the end, stop counting it
            if (pos < x[1]) break;

            // if below interval, update index
            if (pos <= x[pos_index]) {
                pos_index -= 1;
            }

            pos_count(pos_index) += 1;
        }
    }
}


// random numbers with gaussian distribution
double gaussian_deviate(long * idum)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if ( idum < 0) iset =0;
  if (iset == 0) {
    do {
      v1 = 2.*ran2(idum) -1.0;
      v2 = 2.*ran2(idum) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset =0;
    return gset;
  }
} // end function for gaussian deviates

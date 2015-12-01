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

void mc_sampling(int, double, int, double, double, vec&);
// function for gaussian random numbers
double gaussian_deviate(long *);
void output(int, int, vec&, vec&);


int main()
{
    double mc_number = 100000;
    double move_prob=0.5;
    double D=1;
    double T=0.01;
    double d=1;
    double dt=0.00001;
    double l0=sqrt(2*D*dt);
    int Nt = 2; //round(T/dt)+1;
    int Nx = round(d/l0)+1;
    cout << Nt << endl;
    cout << Nx << endl;
    // cout << l0 << endl;
    vec pos_count = zeros<vec>(Nx+1);
    //vec walk2_cumulative = zeros<vec>(Nx+1);

    // Do the mc sampling
    mc_sampling(mc_number, T, Nx, dt, move_prob, pos_count);
    // Print out results
    //output(max_trials, number_walks, walk_cumulative, walk2_cumulative);

    ofstream ofile("/home/vilde/Documents/FYS3150/project5/test2.txt");
    for (int i=0; i<=Nx; i++) {
        ofile << pos_count(i)/(mc_number) << endl;
    }
    ofile.close();

    //cout << "Hello World!" << endl;
    return 0;
}

void mc_sampling(int mc_number, double T, int Nx, double dt, double move_prob, vec& pos_count) {
    long idum = -1;
    pos_count(0)=mc_number; // Number of particles entering the cleft must be constant
    int Nt = 1;
    for (int particles=1; particles<=mc_number; particles++) {
        // Particle's position when it enters the cleft
        int pos_index = 0;

        // Count the particle's position at every time step.
        // This corresponds to storing the movement of particles that had a time interval smaller than T
        // This corresponds to incorporating the number of particles entering the the system at a time later than t=0,
        // having less than time T in the cleft
        for (double time = dt; time <= T; time+=dt) {
        //for (int i=1; i<=Nt; i++) {
            if (ran0(&idum) <= move_prob) { // move to the right

                pos_index += 1;
                if (pos_index >= Nx){        // if it has reached the end, stop observations
                    break;
                }
            } else {                        // move to the left
                pos_index -= 1;
                if (pos_index <=0) {          // if it has returned to start, stop observations
                    break;
                }
            }
            pos_count(pos_index) += 1;
            //walk2_cumulative(walks) += ;
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



void output(int max_trials, int number_walks, vec &walk_cumulative, vec &walk2_cumulative) {
    ofstream ofile("/home/vilde/Documents/FYS3150/project5/testing.txt");
    for( int i = 1; i <= number_walks; i++){
        double xaverage = walk_cumulative(i)/((double) max_trials);
        double x2average = walk2_cumulative(i)/((double) max_trials);
        double variance = x2average - xaverage*xaverage;
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(6) << i;
        ofile << setw(15) << setprecision(8) << xaverage;
        ofile << setw(15) << setprecision(8) << variance << endl;
    }
    ofile.close();
} // end of function output

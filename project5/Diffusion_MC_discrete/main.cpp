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
void output(int, int, vec&, vec&);
void set_up(char**, double&, double&, int&);



/* Main:
 * This program solves the diffusion equation by Monte Carlo simulations and random walk.
 * The position values are discretized.
 * Command line input (in correct order):
 * outfilename, T, dt, mc_number
 */
int main(int argc, char *argv[])
{
    int mc_number = 0;
    double T, dt = 0;
    double move_prob=0.01;
    double D=1;

    set_up(argv, T, dt, mc_number);

    double d=1;
    double l0=sqrt(2*D*dt);
    int Nt = round(T/dt)+1;
    int Nx = round(d/l0)+1;

    vec pos_count = zeros<vec>(Nx+1);

    // Monte Carlo simulation
    mc_sampling(mc_number, T, Nx, dt, move_prob, pos_count);

    // Write to file
    ofstream ofile(argv[1]);
    for (int i=0; i<=Nx; i++) {
        //cout << "hei" << endl;
        ofile << pos_count(i)/(mc_number) << endl;
    }
    ofile.close();

    return 0;
}

void set_up(char** argv, double& T, double& dt, int& mc_number) {
    T = atof(argv[2]);
    dt = atof(argv[3]);
    mc_number = atof(argv[4]);

}


/* Monte Carlo, equal probability:
 * The function solves the diffusion equation by employing Monte Carlo methods and random walks
 * Equal probability of moving a step length l0 to the left or right
 */
void mc_sampling(int mc_number, double T, int Nx, double dt, double move_prob, vec& pos_count) {

    long idum = -1;
    pos_count(0)=mc_number; // Number of particles entering the cleft constant

    for (int particles=1; particles<=mc_number; particles++) {
        // Particle's position when it enters the cleft
        int pos_index = 0;

        /*
         * Count the particle's position at every time step.
         * This corresponds to storing the movement of particles that had a time interval smaller than T
         * This corresponds to incorporating the number of particles entering the system at a time later than t=0,
         * having less than time T in the cleft
         */
        for (double time = dt; time <= T; time+=dt) {
            if (ran0(&idum) <= move_prob) { // move to the right
                pos_index += 1;                
                // if it has reached the end, stop counting it
                if (pos_index >= Nx) break;
            } else {                        // move to the left
                pos_index -= 1;
                // if it has returned to start, stop counting it
                if (pos_index <=0) break;
            }
            // count the particle's position
            pos_count(pos_index) += 1;
        }
    }
}

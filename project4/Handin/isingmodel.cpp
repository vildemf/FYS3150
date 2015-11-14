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

void initialize(int, double, mat&, double&, double&, long&, double, int);
void output(int, int, double, vec&, double, ofstream&, int, vec&, bool);
void metropolis(int, long &, mat &, double &, double &, vec&, int&);
inline int periodic(int, int, int);

/*
 * This program aims at computing the data needed in a) to d).
 * To that end it includes more if-tests making it not so efficient, but capable
 * of producing data with/without thermalization, counting of
 * energy values, ordered/random initial configuration, with a few changes of the
 * initial parameters.
 *
 * For bigger computations on the system, use the parallelized version.
 */

int main(int argc, char* argv[])
{
    // Set up
    ofstream outfile;

    // These values are to be specified by the user to make the program compute the data for a)-d):
    outfile.open("/home/vilde/Documents/FYS3150/project4/b_ordered_L2t1mcE4.txt");
    int L = 2;                    // Number of spins in one direction
    double ord2rand_T = 3.0;      // Alter to force ordered or random config, depending on your T
    int number_mcs = 5E3;         // Number of Monte Carlo cycles
    bool thermalization = false;  // Wait until equilibrium is reached before collecting data?
    bool output_pr_cycle = true;  // Write to file after every MC cycle, rather than one output pr. temperature?
    bool write_Ecount = false;    // This clutters the output file, remove if energy counting not needed
    double start_temp = 0.5;      // Start temperatue
    double end_temp = 1.0;        // End temperature
    double step_temp = 5;         // Step in temperature values


    // set up
    time_t start, finish;     // to measure computing time
    long idum = -1;           // seed for the random numbers
    mat s_mat(L,L);           // spin matrix - one element is one spin

    // Set up temperature array
    int n_temp = round((end_temp-start_temp)/step_temp + 1.0);
    cout << n_temp << endl;
    vec temps = linspace<vec>(start_temp, end_temp, n_temp);
    if ( abs((temps[1]-temps[0])-step_temp) > 1E-5) {
        cout << "Warnng in temperatures - step is a bit different from your input." << endl;
    }

    double temp = 0;      // initialize
    double Cv_prev = 1E5; // used in thermalization estimation
    // for each T
    for (int T_i = 0; T_i < n_temp; T_i++) {
        temp = temps[T_i];
        double E = 0;
        double M = 0;
        vec w = zeros<vec>(17); // 17, to relate deltaE to indices (move out of loop?)
        for( int de =-8; de <= 8; de+=4) {
            w[de+8] = exp(-de/temp);
        }
        vec average = zeros<vec>(5); // 5 is number of relevant expectation values
        vec energies = zeros<vec>(1601);

        initialize(L, temp, s_mat, E, M, idum, ord2rand_T, T_i);
        int a_count = 0; // Counting number of accepted moves
        time(&start);

        for (int cycles = 0; cycles < number_mcs; cycles++){
            metropolis(L, idum, s_mat, E, M, w, a_count);
            if (thermalization) {
                // To check thermalization requirement
                double Erat = (average(0)+E)/average(0);
                double Evariance = (average(1)/(double) number_mcs - average(0)*average(0)/(double)(number_mcs*number_mcs));
                double Cv = Evariance/(temp*temp);
                double Cvrat = Cv/Cv_prev;
                Cv_prev = Cv;
                if (0.95<= Erat<=1.05 && 0.95<=Cvrat<=1.05) {
                    if (L==20) energies(E+800) += 1; // count times each energy appears
                    average(0) += E; average(1) += E*E;
                    average(2) += M; average(3) += M*M; average(4) += fabs(M);
                }
            } else {
                if (L==20) energies(E+800) += 1; // count times each energy appears
                average(0) += E; average(1) += E*E;
                average(2) += M; average(3) += M*M; average(4) += fabs(M);
            }

            if ( output_pr_cycle && (cycles!=number_mcs-1) ) {
                time(&finish);
                output(L, cycles+1, temp, average, difftime(finish, start), outfile, a_count, energies, write_Ecount);
            }
        }

        time(&finish);
        output(L, number_mcs, temp, average, difftime(finish, start), outfile, a_count, energies, write_Ecount);
    } // done with each T

    outfile.close();
    return 0;
}

/* Implements the Metropolis algorithm
 * Proposes L*L=N spins to be flipped
 * New configurations accepted or rejected according to the Metropolis test,
 * if accepted, measurements updated
 */
void metropolis(int L, long &idum, mat &s_mat, double &E, double &M, vec &w, int &a_count) {
    for (int s=0; s<L*L; s++) {
        // Pick random spin for proposed flip
        int x = (int) ( ran1(&idum) * (double)L );
        int y = (int) ( ran1(&idum) * (double)L );
        // Calculate energy change
        int deltaE = 2*s_mat(x,y)*                 // current spin
                (s_mat(x, periodic(y, L, 1)) +     // above spin
                 s_mat(x, periodic(y, L, -1)) +    // below spin
                 s_mat(periodic(x, L, 1), y) +     // to the right spin
                 s_mat(periodic(x, L, -1), y));    // to the left spin
        if ( ran1(&idum) <= w[deltaE + 8] ) {      // if flip accepted-> new configuration
            s_mat(x,y) *= -1;                      // update spin (flip)
            M += (double) 2*s_mat(x,y);              // update M w/ delta M
            E += (double) deltaE;                  // update E w/ delta E
            a_count++;
        }
    }
}



/* Ensures periodic behaviour when reaching end points
 * of the spin matrix
 * inline saves time (at a cost of space)
 * Arguments:
 * - "current" is coordinate of current spin
 * - "L" is number of spins in given direction
 * - "step" is +/-1, indicates which neighbour to move to
 * Returns:
 * - the coordinate of the wanted neighbour spin
 */
inline int periodic(int current, int L, int step) {
    return (current+L+step) % (L);
}




// function to initialise energy, spin matrix and magnetization
void initialize(int L, double temp, mat &s_mat, double &E, double &M, long &idum, double ord2rand_T, int T_i) {
    // setup spin matrix and intial magnetization
    for(int x =0; x <L; x++) {
        for (int y= 0; y <L; y++){
            if (temp < ord2rand_T) {  // Situation 1 - ordered
                s_mat(x, y) = 1; // spin orientation for the ground state
            } else if (T_i==0) {      // Situation 2 - first temp, create random
                if (ran1(&idum) > 0.5) {
                    s_mat(x, y) = 1; // random spin orientation
                } else {
                    s_mat(x, y) = -1;
                }
            }                         // Situation 3 - not first temp, keep random from previous

            M += (double) s_mat(x,y);
        }
    }
    // setup initial energy

    for(int x =0; x < L; x++) {
        for (int y= 0; y < L; y++){
            E -= (double) s_mat(x,y)*
                    (s_mat(x,periodic(y,L,-1)) +
                     s_mat(periodic(x,L,-1), y));
        }
    }
}// end function initialise

void output(int L, int number_mcs, double temp, vec &average, double time, ofstream &outfile, int a_count, vec& energies, bool write_Ecount) {

    // Take MC average
    double Eaverage = average(0)/( (double)number_mcs);
    double E2average = average(1)/( (double)number_mcs);
    double Maverage = average(2)/( (double)number_mcs);
    double M2average = average(3)/( (double)number_mcs);
    double Mabsaverage = average(4)/( (double)number_mcs);

    // Find variances and compute values
    double Evariance = E2average - Eaverage*Eaverage;
    double Mvariance = M2average - Maverage*Maverage;
    double Mabsvariance = M2average - Mabsaverage*Mabsaverage;
    double Cv = Evariance/(temp*temp);
    double X = Mvariance/(temp);
    double X_abs = Mabsvariance/temp;


    // To take values pr number of spins
    double N = (double) L*L;

    outfile << L << "  ";
    outfile << number_mcs << "  ";
    outfile << setprecision(5) << time << "  "; // computation time since start of this temp [sec]
    outfile << setprecision(8) << temp << "  ";
    outfile << setprecision(8) << Eaverage/N << "  ";
    outfile << setprecision(8) << fabs(Maverage)/N << "  ";
    outfile << setprecision(8) << Mabsaverage/N << "  ";
    outfile << setprecision(8) << Cv/N << "  ";
    outfile << setprecision(8) << X/N << "  ";
    outfile << setprecision(8) << X_abs/N << "  ";
    outfile << a_count << "  ";
    outfile << setprecision(8) << Evariance << "  ";
    if (L==20 && write_Ecount) {
        for (int i=0; i<1601; i++) {
            outfile << energies(i) << "  ";
        }
    }
    outfile << endl;
}



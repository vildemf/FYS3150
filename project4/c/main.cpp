#include <iostream>
#include <iomanip>
#include <armadillo>
#include "functions.h"
#include <time.h>
#include "lib.h"

using namespace std;
using namespace arma;

void initialize(int, double, mat&, double&, double&);
void output(int, int, double, vec&, double, ofstream&);
void rand_initialize(int, double, mat&, double&, double&, long&);


int main()
{
    // Set up
    ofstream outfile;
    //outfile.open ("/home/vilde/Documents/FYS3150/project4/c_L20temp24mcE2toE7.txt");
    outfile.open ("/home/vilde/Documents/FYS3150/project4/test2.txt");
    time_t start, finish;     // to measure computing time
    long idum = -1;           // seed for the random numbers
    int L = 20;               // number of spins in one direction
    mat s_mat(L,L);           // spin matrix - one element is one spin


    // Set up MC cycles
    //int start_mcs = 100000;
    //int end_mcs = 1000000;
    //int factor_mcs = 10;

    // Set up temperatures
    double start_temp = 1;
    double end_temp = 2;
    double step_temp = 3;

    // for each T
    for (double temp = start_temp; temp <= end_temp; temp+=step_temp) {
        double E = 0;
        double M = 0;
        vec w = zeros<vec>(17); // 17, to relate deltaE to indices (move out of loop?)
        for( int de =-8; de <= 8; de+=4) {
            w[de+8] = exp(-de/temp);
        }
        vec average = zeros<vec>(5); // 5 is number of relevant expectation values
        //initialize(L, temp, s_mat, E, M);
        rand_initialize(L, temp, s_mat, E, M, idum);


        int cycles = 0;
        time(&start);

        for (cycles = cycles; cycles < 1E1; cycles++){
            metropolis(L, idum, s_mat, E, M, w);
            // update expectation values
            average(0) += E; average(1) += E*E;
            average(2) += M; average(3) += M*M; average(4) += fabs(M);
        }
        time(&finish);
        output(L, cycles, temp, average, difftime(finish, start), outfile);

        for (cycles = cycles; cycles < 1E2; cycles++){
            metropolis(L, idum, s_mat, E, M, w);
            // update expectation values
            average(0) += E; average(1) += E*E;
            average(2) += M; average(3) += M*M; average(4) += fabs(M);
        }
        time(&finish);
        output(L, cycles, temp, average, difftime(finish, start), outfile);


        for (cycles = cycles; cycles < 1E3; cycles++){
            metropolis(L, idum, s_mat, E, M, w);
            // update expectation values
            average(0) += E; average(1) += E*E;
            average(2) += M; average(3) += M*M; average(4) += fabs(M);
        }
        time(&finish);
        output(L, cycles, temp, average, difftime(finish, start), outfile);

        for (cycles = cycles; cycles < 1E4; cycles++){
            metropolis(L, idum, s_mat, E, M, w);
            // update expectation values
            average(0) += E; average(1) += E*E;
            average(2) += M; average(3) += M*M; average(4) += fabs(M);
        }
        time(&finish);
        output(L, cycles, temp, average, difftime(finish, start), outfile);

        for (cycles = cycles; cycles < 1E5; cycles++){
            metropolis(L, idum, s_mat, E, M, w);
            // update expectation values
            average(0) += E; average(1) += E*E;
            average(2) += M; average(3) += M*M; average(4) += fabs(M);
        }
        time(&finish);
        output(L, cycles, temp, average, difftime(finish, start), outfile);

        for (cycles = cycles; cycles < 1E6; cycles++){
            metropolis(L, idum, s_mat, E, M, w);
            // update expectation values
            average(0) += E; average(1) += E*E;
            average(2) += M; average(3) += M*M; average(4) += fabs(M);
        }
        time(&finish);
        output(L, cycles, temp, average, difftime(finish, start), outfile);

    } // done with each T

    outfile.close();

    return 0;
}

// function to initialise energy, spin matrix and magnetization
void initialize(int L, double temp, mat &s_mat, double &E, double &M) {
    // setup spin matrix and intial magnetization
    for(int x =0; x <L; x++) {
        for (int y= 0; y <L; y++){
            if (temp < 4) {
                s_mat(x, y) = 1; // spin orientation for the ground state
            }
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

void rand_initialize(int L, double temp, mat &s_mat, double &E, double &M, long &idum) {
    // setup spin matrix and intial magnetization
    for(int x =0; x <L; x++) {
        for (int y= 0; y <L; y++) {
            if (ran1(&idum) > 0.5) {
                s_mat(x, y) = 1; // random spin orientation
            } else {
                s_mat(x, y) = -1;
            }
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




void output(int L, int number_mcs, double temp, vec &average, double time, ofstream &outfile) {
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
    outfile << setprecision(8) << X_abs/N << endl;

    /*
    cout << "<E> = " << Eaverage/N << endl;
    cout << "|<M>| = " << fabs(Maverage)/N << endl;
    cout << "<|M|> = " << Mabsaverage/N << endl;
    cout << "Cv = " << Cv/N << endl;
    cout << "X = " << X/N << endl;
    cout << "X_abs = " << X_abs/N << endl;
    */
}

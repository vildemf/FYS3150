#include <iostream>
#include <iomanip>
#include <armadillo>
#include "functions.h"
#include <time.h>
#include "lib.h"
#include "mpi.h"

using namespace std;
using namespace arma;

void initialize(int, double, mat&, double&, double&, long&);
void output(int, int, double, vec&, double, ofstream&, int);

int main(int argc, char* argv[])
{
    int number_mcs = 1E5;
    int numprocs, my_rank;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // Set up
    ofstream outfile;
    char filename[10000];
    sprintf(filename, "/home/vilde/Documents/FYS3150/project4/test_results_%d.txt", my_rank);
    outfile.open(filename);
    time_t start, finish;     // to measure computing time
    long idum = -1*my_rank;           // seed for the random numbers
    int L = 60;               // number of spins in one direction
    mat s_mat(L,L);           // spin matrix - one element is one spin


    // Set up temperatures
    double start_temp = 2.0;
    double end_temp = 2.4;
    double step_temp = 0.05;
    int n_temp = round((end_temp-start_temp)/step_temp + 1.0);
    cout << n_temp << endl;
    vec temps = linspace<vec>(start_temp, end_temp, n_temp);
    if ( abs((temps[1]-temps[0])-step_temp) > 1E-10) {
        cout << "Warnng in temps" << endl;
    }
    /*
    for(int i=0; i<numprocs; i++) {
        if(my_rank==i) {
            cout << my_rank << ":" << endl;
            temps.print();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    */
    int my_n = (int) n_temp/numprocs;

    int T_i_start = my_rank*my_n;
    int T_i_stop = (my_rank+1)*my_n;
    if(my_rank == numprocs-1) T_i_stop = n_temp;
    double temp = 0; //initialize

    double Cv_prev = 1E5;
    // for each T
    for (int T_i = T_i_start; T_i < T_i_stop; T_i++) {
        cout << my_rank << " working on T_i=" << T_i << " of " << n_temp << endl;
        temp = temps[T_i];
        double E = 0;
        double M = 0;

        //cout << "rank: " << my_rank << " temp: " << temp << endl;

        vec w = zeros<vec>(17); // 17, to relate deltaE to indices (move out of loop?)
        for( int de =-8; de <= 8; de+=4) {
            w[de+8] = exp(-de/temp);
        }
        vec average = zeros<vec>(5); // 5 is number of relevant expectation values
        vec energies = zeros<vec>(1601);
        //cout << my_rank << ": Initialize..." << endl;
        initialize(L, temp, s_mat, E, M, idum);


        int a_count = 0;
        time(&start);
        //cout << my_rank << ": cycles begin..." << endl;
        for (int cycles = 0; cycles < number_mcs; cycles++){
            //if(my_rank==1) cout << "Cycle: " << cycles << endl;
            metropolis(L, idum, s_mat, E, M, w, a_count, my_rank);

            //if(my_rank==1) cout << "metropolis done." << endl;
            // update expectation values
            //energies(E+800) += 1; // count times each energy appears
            //cout << E << endl;

            // To check thermalization
            double Erat = (average(0)+E)/average(0);
            double Evariance = (average(1)/(double) number_mcs - average(0)*average(0)/(double)(number_mcs*number_mcs));
            double Cv = Evariance/(temp*temp);
            double Cvrat = Cv/Cv_prev;
            Cv_prev = Cv;

            if (0.95<= Erat<=1.05 && 0.95<=Cvrat<=1.05) {
                average(0) += E; average(1) += E*E;
                average(2) += M; average(3) += M*M; average(4) += fabs(M);
            }
        }
        //cout << my_rank << ": cycles done..." << endl;
        time(&finish);
        output(L, number_mcs, temp, average, difftime(finish, start), outfile, a_count);


    } // done with each T

    outfile.close();
    MPI_Finalize ();
    return 0;
}

// function to initialise energy, spin matrix and magnetization
void initialize(int L, double temp, mat &s_mat, double &E, double &M, long &idum) {
    // setup spin matrix and intial magnetization
    for(int x =0; x <L; x++) {
        for (int y= 0; y <L; y++){
            if (temp < 1.5) {
                s_mat(x, y) = 1; // spin orientation for the ground state
            } else {
                if (ran1(&idum) > 0.5) {
                    s_mat(x, y) = 1; // random spin orientation
                } else {
                    s_mat(x, y) = -1;
                }
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
/*
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
*/
void output(int L, int number_mcs, double temp, vec &average, double time, ofstream &outfile, int a_count) {
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
    outfile << setprecision(8) << Evariance << endl;

    /*
    cout << "<E> = " << Eaverage/N << endl;
    cout << "|<M>| = " << fabs(Maverage)/N << endl;
    cout << "<|M|> = " << Mabsaverage/N << endl;
    cout << "Cv = " << Cv/N << endl;
    cout << "X = " << X/N << endl;
    cout << "X_abs = " << X_abs/N << endl;
    */
}

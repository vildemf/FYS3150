#include <iostream>
#include <armadillo>
#include "functions.h"
#include <time.h>

using namespace std;
using namespace arma;

void initialize(int, double, mat&, double&, double&);

int main()
{
    // before T
    long idum = -1;
    int L = 2;
    int number_mcs = 1000000;
    mat s_mat(L,L);
    //int cc = 0;
    // for each T
    double temp = 1.0;
    double E = 0;
    double M = 0;
    vec w = zeros<vec>(17); // 17 is convenient in order to relate deltaE to indices
    for( int de =-8; de <= 8; de+=4) {
        w[de+8] = exp(-de/temp);
     }
    vec average = zeros<vec>(5); // 5 is number of relevant expectation values
    initialize(L, temp, s_mat, E, M);
    for (int cycles = 1; cycles <= number_mcs; cycles++){
        metropolis(L, idum, s_mat, E, M, w);        
        //if (abs(E+8)>1e-3) {
        //    cout << E << " " << ++cc << endl;
        //}
        // update expectation values
        average(0) += E; average(1) += E*E;
        average(2) += M; average(3) += M*M; average(4) += fabs(M);
    }
    // take MC average
    average *= 1.0/((double) number_mcs);
    double Eaverage = average(0);
    double E2average = average(1);
    double Maverage = average(2);
    double M2average = average(3);
    double Mabsaverage = average(4);

    // find variances and compute values
    double Evariance = E2average - Eaverage*Eaverage;
    double Mvariance = M2average - Maverage*Maverage;
    double Mabsvariance = M2average - Mabsaverage*Mabsaverage;
    double Cv = Evariance/(temp*temp);
    double X = Mvariance/(temp);
    double X_abs = Mabsvariance/temp;


    // To take values pr number of spins
    double N = (double) L*L;

    //cout << "<M2> = " << M2average << endl;
    //cout << "<E2> = " << E2average << endl;
    cout << "<E> = " << Eaverage/N << endl;
    cout << "|<M>| = " << fabs(Maverage)/N << endl;
    cout << "<|M|> = " << Mabsaverage/N << endl;
    cout << "Cv = " << Cv/N << endl;
    cout << "X = " << X/N << endl;
    cout << "X_abs = " << X_abs/N << endl;

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


//void output(int L, int number_mcs, double temp, vec &average) {
//
//}

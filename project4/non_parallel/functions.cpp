#include <armadillo>
#include <cmath>
#include "lib.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include "functions.h"

using namespace arma;


/* Implements the Metropolis algorithm
 * Proposes L*L=N spins to be flipped
 * New configurations accepted or rejected,
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

/*
void metropolis(int L, long &idum, mat &s_mat, double &E, double &M, vec &w) {
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
        }
    }
}
*/

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
//int periodic(int current, int L, int step) {
//    return (current+L+step) % (L);
//}



#ifndef FUNCTIONS
#define FUNCTIONS

#include <armadillo>
using namespace arma;

void metropolis(int, long &, mat &, double &, double &, vec &, int&, int my_rank);
inline int periodic(int current, int L, int step) {
    return (current+L+step) % (L);
}

#endif // FUNCTIONS


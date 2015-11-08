#ifndef FUNCTIONS
#define FUNCTIONS

#include <armadillo>
using namespace arma;

void metropolis(int, long &, mat &, double &, double &, vec &, int&);
inline int periodic(int current, int L, int step);


#endif // FUNCTIONS


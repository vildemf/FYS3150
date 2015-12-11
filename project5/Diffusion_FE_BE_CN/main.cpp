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


void crank_nicolson(int, int, double, vec&, vec&, vec&, ofstream&, bool);
void backward_euler(int, int, double, vec&, vec&, vec&, ofstream&, bool);
void forward_euler(int, int, double, vec&, vec&, vec&, ofstream&, bool);
double initial(double);
void set_up(double&, double&, double&, double&, double&,
            double&, int&, int&, vec&, vec&, vec&,
            vec&, char**, bool&, bool&, bool&, bool&);


/* Main:
 * Function to solve the diffusion equation using either of three methods: explicit Forward Euler,
 * implicit Backward Euler or implicit Crank-Nicolson.
 *
 * Input on command line in correct order (0/1 indicates false/true):
 * filename, T, d, Nt, Nx, forward=0or1, backward=0or1, crank=0or1, store_x_not_t=0or1
 */
int main(int argc, char *argv[])
{
    // Declerations
    ofstream outfile;
    outfile.open(argv[1]);
    double T, d, D, dt, dx, alpha = 0;
    int Nt, Nx= 0;
    vec t, x, v_prev, v;
    bool forward, backward, crank, store_x_not_t;

    // Initialize
    set_up(T, d, D, dt, dx, alpha, Nt, Nx, t, x, v_prev, v, argv,
           forward, backward, crank, store_x_not_t);

    // Solve with method of choice
    if (forward)  forward_euler(Nt, Nx, alpha, v_prev, v, x, outfile, store_x_not_t);
    if (backward) backward_euler(Nt, Nx, alpha, v_prev, v, x, outfile, store_x_not_t);
    if (crank)    crank_nicolson(Nt, Nx, alpha, v_prev, v, x, outfile, store_x_not_t);

    // Write to file
    if (store_x_not_t) {
        for (int i=0; i<=Nx; i++) {
            outfile << setprecision(8) << v(i) << endl;
        }
    }
    outfile.close();

    return 0;
}



/* Input and set up:
 * The function reads command line input and
 * initializes the variables
 */
void set_up(double& T, double& d, double& D, double& dt, double& dx,
            double& alpha, int& Nt, int& Nx, vec& t, vec& x, vec& v_prev,
            vec& v, char** argv, bool& forward, bool& backward, bool& crank, bool& store_x_not_t) {
    T = atof(argv[2]);        // t(Nt)
    d = atof(argv[3]);        // x(Nx)
    Nt = atof(argv[4]);       // max index for t-values
    Nx = atof(argv[5]);       // max index for x-values

    // Set the choices to true/false values: 1=true, 0=false
    // Choice of method:
    forward = *argv[6]!='0';
    backward = *argv[7]!='0';
    crank = *argv[8]!='0';

    // true = output for x-values, v(x, t=T)
    // false = output for t-values, v(x=0.5, t)
    store_x_not_t = *argv[9]!='0';

    t = linspace<vec>(0, T, Nt+1);
    x = linspace<vec>(0, d, Nx+1);
    dt = t(1) - t(0);
    dx = x(1) - x(0);

    D = 1; // diffusuion coefficient
    alpha = D*dt/(dx*dx);

    v_prev = zeros<vec>(Nx+1); // store previous for computing for next t
    v = zeros<vec>(Nx+1);      // output solution
}


/* Implicit Crank-Nicolson method:
 * Function solves the diffusion equation with implicit Crank-Nicolson method
 * as a linear equation problem
 */
void crank_nicolson(int Nt, int Nx, double alpha, vec&v_prev, vec& v, vec& x, ofstream& outfile, bool store_x_not_t) {
    /* Initialize:
     * Set up the known constants in the matrix
     * and arrays to store row operations
     * for the Backward Euler part of the procedure
     */
    double b = 2+2*alpha; // diagonal
    double a_c = -alpha;  // not diagonal
    double ac = a_c*a_c;  // useful precalculation
    vec v_prev_rowoperated = zeros<vec>(Nx+1); // right side after row operation
    vec b_rowoperated = zeros<vec>(Nx+1);      // diagonal after row operation

    /* Initial:
     * Set initial condition v(x, t=0)
     * (I keep the boundary conditions, v(0),v(Nx)=0, not =initial(), at t=0)
     */
    for (int i=1; i<Nx; i++) {
        v_prev(i) = initial(x(i));
    }

    if (!store_x_not_t) outfile << 0 << " " << v_prev(100) << " " << x(100) << endl;

    /* Computations:
     * For each t, solve for v combining the Forward and Backward Euler procedures.
     * The matrix system excludes the known v(0,t), v(d,t)=0 i=0, i=Nx.
     * Thus we initialize v(xd, dt) i=1, j=1, and work with i=2,...,Nx-1
     */
    for (int j=1; j<=Nt; j++) {

        /* Step 1 (Forward Euler procedure)
         * Solve:
         * V_prev_tilde=(2I-alphaB)V_prev
         */
        for (int i=1; i<Nx; i++) {
            v(i) = alpha*v_prev(i-1) + (2-2*alpha)*v_prev(i) + alpha*v_prev(i+1);
        }
        for (int i=0; i<=Nx; i++) {
            v_prev(i) = v(i);
        }


        /* Step 2 (Backward Euler procedure)
         * Solve:
         * (2I+alphaB)V=V_prev_tilde
         * The matrix system excludes the known v(0,t), v(d,t)=0 i=0, i=Nx.
         * Thus we initialize v(xd, dt) i=1, j=1, and work with i=2,...,Nx-1
         */

        // Forward substitution
        v_prev_rowoperated(1) = v_prev(1);
        b_rowoperated(1) = b;
        for (int i=2; i<Nx; i++) {
            b_rowoperated(i) = b - ac/b_rowoperated(i-1);
            v_prev_rowoperated(i) = v_prev(i) - a_c*v_prev_rowoperated(i-1)/b_rowoperated(i-1);
        }

        // Backward substitution
        v(Nx-1) = v_prev_rowoperated(Nx-1)/b_rowoperated(Nx-1);
        for (int i=Nx-2; i>0; i--) {
            v(i) = (v_prev_rowoperated(i)-a_c*v(i+1))/b_rowoperated(i);
        }

        if (!store_x_not_t) outfile << j << " " << v(100) << endl;

        // Update for new j
        for (int i=0; i<=Nx; i++) {
            v_prev(i) = v(i);
        }
    }
}



/* Implicit Backward Euler method:
 * Function solves the diffusion equation with implicit Backward Euler
 * as a linear equation problem
 */
void backward_euler(int Nt, int Nx, double alpha, vec&v_prev, vec& v, vec& x, ofstream& outfile, bool store_x_not_t) {
    /* Initialize:
     * Set up the known constants in the matrix
     * and arrays to store row operations
     */
    double b = 1+2*alpha; // diagonal
    double a_c = -alpha;  // non diagonals
    double ac = a_c*a_c;  // useful precalculation
    vec v_prev_rowoperated = zeros<vec>(Nx+1); // right side after row operation
    vec b_rowoperated = zeros<vec>(Nx+1);      // diagonal after row operation

    /* Initial:
     * Set initial condition v(x, t=0)
     * (I keep the boundary conditions, v(0),v(Nx)=0, not =initial(), at t=0)
     */
    for (int i=1; i<Nx; i++) {
        v_prev(i) = initial(x(i));
    }

    if (!store_x_not_t) outfile << 0 << " " << v_prev(100) << " " << x(100) << endl;

    /* Computations:
     * For each t, solve the matrix system to find v - forward substitution and backward substitution.
     * The matrix system excludes the known v(0,t), v(d,t)=0 i=0, i=Nx.
     * Thus we initialize v(xd, dt) i=1, j=1, and work with i=2,...,Nx-1
     */
    for (int j=1; j<=Nt; j++) {

        // Forward substitution
        v_prev_rowoperated(1) = v_prev(1);
        b_rowoperated(1) = b;
        for (int i=2; i<Nx; i++) {
            b_rowoperated(i) = b - ac/b_rowoperated(i-1);
            v_prev_rowoperated(i) = v_prev(i) - a_c*v_prev_rowoperated(i-1)/b_rowoperated(i-1);
        }

        // Backward substitution
        v(Nx-1) = v_prev_rowoperated(Nx-1)/b_rowoperated(Nx-1);
        for (int i=Nx-2; i>0; i--) {
            v(i) = (v_prev_rowoperated(i)-a_c*v(i+1))/b_rowoperated(i);
        }

        if (!store_x_not_t) outfile << j << " " << v(100) << endl;

        // Update for new j
        for (int i=0; i<=Nx; i++) {
            v_prev(i) = v(i);
        }
    }

}



/* Explicit Forward Euler method:
 * Function solves the diffusion equation employing explicit Forward Euler
 *
 * Stability condition Forward Euler: dt/(dx**2) < 0.5
 * corresponding to Nt > 2T(Nx**2)
 */
void forward_euler(int Nt, int Nx, double alpha, vec& v_prev, vec& v, vec& x, ofstream& outfile, bool store_x_not_t) {
    /* Initial:
     * Set tinitial condition v(x, t=0)
     * (I keep the boundary conditions, v(0),v(Nx)=0, not =initial(), at t=0)
     */
    for (int i=1; i<Nx; i++) {
        v_prev(i) = initial(x(i));
    }

    if (!store_x_not_t) outfile << 0 << " " << v_prev(100) << " " << x(100) << endl;

    /* Computations:
     * Calculate v(x, t) for t=dt,...,T and x=dx,..,d-dx
     * Corresponds to j=1,...,Nt and i=1,...,Nx-1
     * Leave out i=0 and i=Nx to keep v here 0 (boundary condition)
     */
    for (int j=1; j<=Nt; j++) {
        for (int i=1; i<Nx; i++) {
            v(i) = alpha*v_prev(i-1) + (1-2*alpha)*v_prev(i) + alpha*v_prev(i+1);
        }       
        for (int i=0; i<=Nx; i++) {
            v_prev(i) = v(i);
        }
        if (!store_x_not_t) outfile << j << " " << v(100) << endl;

    }
    if (!store_x_not_t) outfile.close();
}



/*
 * Sets  the initial condition v(x, 0)
 */
double initial(double x) {
    return x-1;
}

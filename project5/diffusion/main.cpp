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


void crank_nicolson(int, int, double, vec&, vec&, vec&);
void backward_euler(int, int, double, vec&, vec&, vec&);
void forward_euler(int, int, double, vec&, vec&, vec&);
double initial(double);

int main(int argc, char *argv[])
{
    if (argc < 3) {
        printf("LOL\n");
        printf("%s", argv[1]);
        exit(EXIT_SUCCESS);
    }
    ofstream outfile;
    outfile.open("/home/vilde/Documents/FYS3150/project5/crank.txt");
    double T = 0.5; // t(Nt)
    double d = 1;  // x(Nx)
    int Nt = 50000;   // number of t-values
    int Nx = 200;   // number of x-values
    vec t = linspace<vec>(0, T, Nt+1);
    vec x = linspace<vec>(0, d, Nx+1);
    double dt = t(1) - t(0);
    double dx = x(1) - x(0);

    double D = 1;
    double alpha = D*dt/(dx*dx);

    vec v_prev = zeros<vec>(Nx+1); // store previous for computing for next t
    vec v = zeros<vec>(Nx+1);

    forward_euler(Nt, Nx, alpha, v_prev, v, x);
    //backward_euler(Nt, Nx, alpha, v_prev, v, x);
    crank_nicolson(Nt, Nx, alpha, v_prev, v, x);

    outfile << Nt << " " << Nx << " " << dt << " " << dx << " " << T << endl;
    for (int i=0; i<=Nx; i++) {
        outfile << setprecision(8) << v(i) << endl;
    }
    outfile.close();

    cout << "Hello World!" << endl;
    return 0;
}


void crank_nicolson(int Nt, int Nx, double alpha, vec&v_prev, vec& v, vec& x) {
    double b = 2+2*alpha; // diagonal
    double a_c = -alpha; // not diagonal
    double ac = a_c*a_c;

    vec v_prev_rowoperated = zeros<vec>(Nx+1);
    vec b_rowoperated = zeros<vec>(Nx+1); // new diagonal b_rowoperated

    // v(x, t=0) j=0
    for (int i=1; i<Nx; i++) {
        v_prev(i) = initial(x(i));
    }

    for (int j=1; j<=Nt; j++) {
        // Step 1: Solve V_prev_tilde = (2I-alphaB)V_prev (Forward Euler procedure)
        for (int i=1; i<Nx; i++) {
            v(i) = alpha*v_prev(i-1) + (2-2*alpha)*v_prev(i) + alpha*v_prev(i+1);
        }
        for (int i=0; i<=Nx; i++) {
            v_prev(i) = v(i);
        }

        // Step 2: Solve (2I+alphaB)V=V_prev_tilde

        // Rember the matrix system excludes the known v(0,t), v(d,t)=0 i=0, i=Nx.
        // Thus we initialize v(xd, dt) i=1, j=1, and work with i=2,...,Nx-1
        v_prev_rowoperated(1) = v_prev(1);
        b_rowoperated(1) = b;
        // Forward sub
        for (int i=2; i<Nx; i++) {
            b_rowoperated(i) = b - ac/b_rowoperated(i-1);
            v_prev_rowoperated(i) = v_prev(i) - a_c*v_prev_rowoperated(i-1)/b_rowoperated(i-1);
        }

        // backward substitution
        v(Nx-1) = v_prev_rowoperated(Nx-1)/b_rowoperated(Nx-1);
        for (int i=Nx-2; i>0; i--) {
            v(i) = (v_prev_rowoperated(i)-a_c*v(i+1))/b_rowoperated(i);
        }

        // update for new j
        for (int i=0; i<=Nx; i++) {
            v_prev(i) = v(i);
        }
    }

}


void backward_euler(int Nt, int Nx, double alpha, vec&v_prev, vec& v, vec& x) {
    double b = 1+2*alpha; // diagonal
    double a_c = -alpha; // not diagonal
    double ac = a_c*a_c;

    vec v_prev_rowoperated = zeros<vec>(Nx+1);
    vec b_rowoperated = zeros<vec>(Nx+1); // new diagonal b_rowoperated

    // v(x, t=0) j=0
    for (int i=1; i<Nx; i++) {
        v_prev(i) = initial(x(i));
    }
    for (int j=1; j<=Nt; j++) {
        // Rember the matrix system excludes the known v(0,t), v(d,t)=0 i=0, i=Nx.
        // Thus we initialize v(xd, dt) i=1, j=1, and work with i=2,...,Nx-1

        v_prev_rowoperated(1) = v_prev(1);
        b_rowoperated(1) = b;

        // Forward sub
        for (int i=2; i<Nx; i++) {
            b_rowoperated(i) = b - ac/b_rowoperated(i-1);
            v_prev_rowoperated(i) = v_prev(i) - a_c*v_prev_rowoperated(i-1)/b_rowoperated(i-1);
        }

        // backward substitution
        v(Nx-1) = v_prev_rowoperated(Nx-1)/b_rowoperated(Nx-1);
        for (int i=Nx-2; i>0; i--) {
            v(i) = (v_prev_rowoperated(i)-a_c*v(i+1))/b_rowoperated(i);
        }

        // update for new j
        for (int i=0; i<=Nx; i++) {
            v_prev(i) = v(i);
        }
    }

}

/*
 * Stability condition forward_euler: dt/(dx**2) < 0.5
 * gives the condition Nt > 2T(Nx**2)
 *
 */
void forward_euler(int Nt, int Nx, double alpha, vec& v_prev, vec& v, vec& x) {
    //ofstream outfile;
    //outfile.open("/home/vilde/Documents/FYS3150/project5/reler.txt");
    // t=0
    // v(0), v(Nx) = 0 or initial() at t=0?
    for (int i=1; i<Nx; i++) {
        v_prev(i) = initial(x(i));
    }

    //outfile << 0 << " " << v_prev(100) << " " << x(100) << endl;
    // Calculate v(x, t) for t=dt,...,T and x=dx,..,d-dx
    // Corresponds to j=1,...,Nt and i=1,...,Nx-1
    // Leave out i=0 and i=Nx to keep v here 0 (boundary condition)
    for (int j=1; j<=Nt; j++) {
        for (int i=1; i<Nx; i++) {
            v(i) = alpha*v_prev(i-1) + (1-2*alpha)*v_prev(i) + alpha*v_prev(i+1);
        }       
        for (int i=0; i<=Nx; i++) {
            v_prev(i) = v(i);
        }
        //outfile << j << " " << v(100) << endl;

    }
    //outfile.close();
}


double initial(double x) {
    return x-1;
}

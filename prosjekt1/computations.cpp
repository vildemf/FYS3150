#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

int main()
{
    int n = 100000;
    double x_0 = 0;
    double x_n_1 = 1;
    double h = (x_n_1 - x_0)/(n+1); // step
    double h2 = h*h;

    double v_0 = 0;
    double v_n_1 = 0;

    double b = 2;                      // diagonal
    double a_c = -1;                   // not diagonal
    double ac = a_c*a_c;
    
    double x[n+2]; // n+2 elements, from 0 to n+1
    for (int i=0; i<n+2; i++)
    {
      x[i] = i*h;
    }
    
    // the linear equation system has elements i=1,...,n (i=0,...,n-1 in array indexing)
    double f[n]; // to store f multiplied by h**2 (so "f" is not an optimal name choice)
    double f_new[n];
    double b_new[n]; // new diagonal
    double v[n]; // exclude v0 and vn+1 from the system since they're zero

    b_new[0] = b;
    f_new[0] = f[0] = h2*100*exp(-10*x[1]); // remember the system excludes f(x0), so f0=f(x1)

    // Forward sub
    int xi; // x needs a seperate index since x0 is included in the array but excluded in the system
    for (int i=1; i<n; i++)
    {
        b_new[i] = b - ac/b_new[i-1];

        xi = i + 1;
        f[i] = h2*100*exp(-10*x[xi]); // f multiplied by h**2

        f_new[i] = f[i] - a_c*f_new[i-1]/b_new[i-1];
    }

    v[n-1] = f_new[n-1]/b_new[n-1];
    for (int i=n-2; i>=0; i--)
    {
        v[i] = (f_new[i]-a_c*v[i+1])/b_new[i];
    }

    ofstream outfile;
    outfile.open ("n10_5.txt");
    for (int i=0; i<n; i++)
    {
      outfile << v[i] << "\n";

      //cout << x[i] << "\n";
      //cout << f_new[i] << "\n";
      //cout << v[i] << "\n";
    }
    outfile.close();    
    
    return 0;
}

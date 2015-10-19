#include <iostream>
#include <armadillo>
#include "legendre.cpp"
#include <time.h>

using namespace std;
using namespace arma;

int main()
{
    // General
    double a = -3.0;
    double b = 3.0;
    time_t start, finish;

    /*
    int N = 20;
    // Legendre
    time(&start);
    double int_gauss = legendre(N, a, b);
    time(&finish);
    cout << "Integral Legendre: " << int_gauss << " Time: "<< difftime(finish, start) << endl;

    // Legendre and Laguerre
    time(&start);
    double int_legendre_laguerre = legendre_and_laguerre(N);
    time(&finish);
    cout << "Integral Legendre and Laguerre: " << int_legendre_laguerre << " Time: " << difftime(finish, start) << endl;
    */

    int n = 1.0e8;
    // Brute force Monte Carlo
    a = -5;
    b = 5;
    double mu = 0;
    double variance = 0;
    /*
    time(&start);
    brute_MC(n, a, b, mu, variance);
    time(&finish);
    double std = sqrt(variance);
    cout << "Integral brute force MC: " << mu << " Standard deviation: " << std << " Time: " << difftime(finish, start) << endl;
    */

    // Improved Monte Carlo using importance sampling
    double u_b = 10;
    mu = 0;
    variance = 0;
    time(&start);
    improved_MC(n, u_b, mu, variance);
    time(&finish);
    double std = sqrt(variance);
    cout << "Integral brute force MC: " << mu << " Standard deviation: " << std << " Time: " << difftime(finish, start) << endl;



    return 0;
}





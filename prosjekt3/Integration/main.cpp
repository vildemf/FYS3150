#include <iostream>
#include <armadillo>
#include "methods.cpp"
#include <time.h>

using namespace std;
using namespace arma;

int main()
{
    // General
    double a = -3.0;
    double b = 3.0;
    time_t start, finish;
    double exact = 5*M_PI*M_PI/(16*16);
    ofstream outfile0;
    outfile0.open ("/home/vilde/Documents/FYS3150/prosjekt3/output.txt");
    /*
    int N = 35;
    // Legendre
    time(&start);
    double int_gauss = legendre(N, a, b);
    time(&finish);
    cout << "N= " << N << " Integral Legendre: " << int_gauss <<
            " Relative error: " << abs(exact-int_gauss)/exact <<
            " Time: "<< difftime(finish, start) << endl;
    outfile0 << "N= " << N << " Integral Legendre: " << int_gauss <<
            " Relative error: " << abs(exact-int_gauss)/exact <<
            " Time: "<< difftime(finish, start) << endl;


    // Legendre and Laguerre
    time(&start);
    double int_legendre_laguerre = legendre_and_laguerre(N);
    time(&finish);
    cout << "N= " << N << " Integral Legendre and Laguerre: " << int_legendre_laguerre <<
            " Relative error: " << abs(exact-int_legendre_laguerre)/exact << " Time: " << difftime(finish, start) << endl;
    outfile0 << "N= " << N << " Integral Legendre and Laguerre: " << int_legendre_laguerre <<
            " Relative error: " << abs(exact-int_legendre_laguerre)/exact << " Time: " << difftime(finish, start) << endl;
    */

    int n = 1.0e8;
    // Brute force Monte Carlo
    a = -3;
    b = 3;
    double mu = 0;
    double variance = 0;

    time(&start);
    brute_MC(n, a, b, mu, variance);
    time(&finish);
    double std = sqrt(variance)/sqrt(n);
    cout << "N= " << n << " Integral brute force MC: "
        << mu << " Standard deviation: " << std <<
           " Relative error: " << abs(exact-mu)/exact <<
           " Time: " << difftime(finish, start) << endl;
    outfile0 << "N= " << n << " Integral brute force MC: "
        << mu << " Standard deviation: " << std <<
           " Relative error: " << abs(exact-mu)/exact <<
           " Time: " << difftime(finish, start) << endl;



    // Improved Monte Carlo using importance sampling

    double u_b = 10;
    mu = 0;
    variance = 0;
    time(&start);
    improved_MC(n, u_b, mu, variance);
    time(&finish);
    std = sqrt(variance)/sqrt(n);
    cout << "N= "<< n <<" Integral improved MC: "
            << mu << " Standard deviation: " << std <<
               " Relative error: " << abs(exact-mu)/exact <<
               " Time: " << difftime(finish, start) << endl;

    outfile0 << "N= "<< n <<" Integral improved MC: "
            << mu << " Standard deviation: " << std <<
               " Relative error: " << abs(exact-mu)/exact <<
               " Time: " << difftime(finish, start) << endl;


    outfile0.close();
    return 0;
}





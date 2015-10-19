#include <armadillo>
#include <cmath>
//#include <math.h>
#include "lib.h"

#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10
//#define PI 3.14159265359



using namespace arma;

void improved_MC(int n, double u_b, double &mu, double &variance);
double int_func_MC(double theta1, double theta2, double phi1, double phi2, double u1, double u2);
void brute_MC(int n, double a, double b, double &mu, double &var);
double int_function1(double xi, double xj, double xk, double xl, double xm, double xn);
double int_function2(double theta1, double theta2, double phi1, double phi2, double u1, double u2);
double legendre(int N, double a, double b);
void gauss_laguerre(double *x, double *w, int n, double alf);
double gammln( double xx);
void gauleg1(double a, double b, double x[], double w[], int N);



void improved_MC(int n, double u_b, double &mu, double &variance) {
    long idum = time(0);
    double sum_sigma, f;
    sum_sigma =f = 0.;
    double dimension = 6;
    vec x(dimension);
    double theta1, theta2, phi1, phi2, u1, u2;
    theta1=theta2=phi1=phi2=u1=u2 = 0;
    double a_theta, b_theta, a_phi, b_phi, a_u;
    a_theta = a_phi = a_u = 0;
    b_theta = M_PI;
    b_phi = 2*M_PI;

    double alpha = 2.0;
    double jacdet_theta = pow(b_theta - a_theta, 2);
    double jacdet_phi = pow(b_phi - a_phi, 2);

    for (int i=0; i<n; i++) {

        for (int j=0; j<dimension; j++) {
            x(j)=ran0(&idum);
        }
        // Uniform dist. mapping
        theta1 = x(0)*(b_theta-a_theta) + a_theta;
        theta2 = x(1)*(b_theta-a_theta) + a_theta;
        phi1 = x(2)*(b_phi-a_phi) + a_phi;
        phi2 = x(3)*(b_phi-a_phi) + a_phi;
        // Exp. dist. mapping
        u1 = -log(1-x(4));
        u2 = -log(1-x(5));
        f = int_func_MC(theta1, theta2, phi1, phi2, u1, u2);
        mu += f;
        sum_sigma += f*f;
    }
    double factor = pow(2*alpha, -5)*jacdet_phi*jacdet_theta;
    mu = factor*mu/((double) n );
    sum_sigma = pow(factor, 2)*sum_sigma/((double) n );
    variance=sum_sigma-mu*mu;
}

// To be used by the improved MC
double int_func_MC(double theta1, double theta2, double phi1, double phi2, double u1, double u2) {
    double cosBeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    double u12_sq = u1*u1 + u2*u2 - 2*u1*u2*cosBeta;
    if (abs(u12_sq)<pow(10., -6) || u12_sq<0) { // Can never happen physically, so skip
        return 0;
    }
    // exp(-(u1+u2)) is accounted for by dx/du, here is F(y(x)):
    double function_value = sin(theta1)*sin(theta2)*u1*u1*u2*u2/sqrt(u12_sq);
    return function_value;
}


void brute_MC(int n, double a, double b, double &mu, double &variance) {
    long idum = time(0);
    double x, sum_sigma, f;
    x =sum_sigma =f =0.;
    int dimension = 6; // dimension of the integral
    double jacdet = pow((b-a), dimension);
    //int j = 0;
    //vec y(dimension);
    double *y = new double[dimension];
    // evaluate the integral with the a crude Monte-Carlo method
    for (int i = 0; i < n; i++){
        //cout << "i = " << i << endl;
        for (int j = 0; j<dimension; j++) {// generate 6 random points to evaluate
            x=ran0(&idum); // random number with limits [0, 1]
            //cout << "x = " << x;
            y[j] = a + (b-a)*x; // change of variables to fit our limits, by transformed uniform distribution
            //cout << "  y = " << y[j] << endl;
        }
        f=int_function1(y[0], y[1], y[2], y[3], y[4], y[5]);
        mu += f;
        sum_sigma += f*f;
    }
    mu = jacdet*mu/((double) n );
    sum_sigma = pow(jacdet, 2)*sum_sigma/((double) n );
    variance=sum_sigma-mu*mu;

}

double int_function1(double xi, double xj, double xk, double xl, double xm, double xn) {
    double alpha = 2.0; // Charge of the Helium atom
    double relative_distance_sq = pow((xi-xl), 2) + pow((xj-xm), 2) + pow((xk-xn), 2);
    double r1 = sqrt(xi*xi + xj*xj + xk*xk);
    double r2 = sqrt(xl*xl + xm*xm + xn*xn);

    if (abs(relative_distance_sq)<pow(10., -6) || relative_distance_sq<0) { // Can never happen physically, so skip
        return 0;
    }
    double function_value = (1.0/sqrt(relative_distance_sq))*exp(-2*alpha*(r1 + r2));
    return function_value;
}

double int_function2(double theta1, double theta2, double phi1, double phi2, double u1, double u2) {
    // Need to multiply the final integral by (2*alpha)**(-5)!! It's not here
    //double alpha = 2.0;
    double cosBeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    double u12_sq = u1*u1 + u2*u2 - 2*u1*u2*cosBeta;
    if (abs(u12_sq)<pow(10., -6) || u12_sq<0) { // Can never happen physically, so skip
        return 0;
    }
    // exp(-(u1+u2))*pow(u1,alf)*pow(u2,alf) is accounted for by W, here is g:
    double function_value = sin(theta1)*sin(theta2)/sqrt(u12_sq);
    return function_value;
}


double legendre_and_laguerre(int N) {
    double alf = 2.0;
    double alpha = 2.0;
    double *theta = new double[N];
    double *wtheta = new double[N];
    double *phi = new double[N];
    double *wphi = new double[N];
    double *u1 = new double[N+1];
    double *wu1 = new double[N+1];
    gauleg(0, 2*M_PI, phi, wphi, N);
    gauleg(0, M_PI, theta, wtheta, N);
    gauss_laguerre(u1, wu1, N, alf);
    double int_legendre_laguerre = 0;
    for (int i=0;i<N;i++){
       for (int j = 0;j<N;j++){
           for (int k = 0;k<N;k++){
               for (int l = 0;l<N;l++){
                   for (int m = 1;m<N+1;m++){
                       for (int n = 1;n<N+1;n++){
                           int_legendre_laguerre+=wtheta[i]*wtheta[j]*wphi[k]*wphi[l]*wu1[m]*wu1[n]
                                   *int_function2(theta[i],theta[j],phi[k],phi[l],u1[m],u1[n]);
                       }
                   }
               }
           }
       }

    }
    return int_legendre_laguerre/pow(2*alpha, 5); // included the alpha factor excluded in the integrand
}

double legendre(int N, double a, double b) // integrand func
{
// array for integration points and weights using Legendre polynomials
     double *x = new double [N];
     double *w = new double [N];
//   set up the mesh points and weights
     //cout << "enter gauleg" << endl;
     gauleg(a,b,x,w, N);
//   evaluate the integral with the Gauss-Legendre method
//   Note that we initialize the sum
     double int_legendre = 0.;
//   six-double loops
     //for (int i = 0; i<N; i++) {
     //    cout << "i = " << i << " x = " << x[i] << " w = " << w[i] << endl;
     //}

     //cout << "Enter integration loop" << endl;
     for (int i=0;i<N;i++){
        for (int j = 0;j<N;j++){
            for (int k = 0;k<N;k++){
                for (int l = 0;l<N;l++){
                    for (int m = 0;m<N;m++){
                        for (int n = 0;n<N;n++){
                            int_legendre+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function1(x[i],x[j],x[k],x[l],x[m],x[n]);
                        }
                    }
                }
            }
        }
     }
     return int_legendre;
}



void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;
    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
            (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}
// end function gaulag

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}



void gauleg1(double x1, double x2, double x[], double w[], int n)
{
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;
    double const pi = 3.14159265359;
    double *x_low, *x_high, *w_low, *w_high;
    m = (n + 1)/2; // roots are symmetric in the interval
    xm = 0.5 * (x2 + x1);
    //cout << "m: " << m << endl;
    xl = 0.5 * (x2 - x1);
    x_low = x; // pointer initialization
    x_high = x + n - 1;
    w_low = w;
    w_high = w + n - 1;
    //cout << "Enter for loop" << endl;
    for(i = 1; i <= m; i++) { // loops over desired roots
        z = cos(pi * (i - 0.25)/(n + 0.5));
        /*
        ** Starting with the above approximation to the ith root
        ** we enter the mani loop of refinement bt Newtons method.
        */
        //cout << "New do while, i = " << i << " of " << m << endl;
        do {
            p1 =1.0;
            p2 =0.0;
            /*
            ** loop up recurrence relation to get the
            ** Legendre polynomial evaluated at x
            */
            for(j = 1; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
            }
            /*
            ** p1 is now the desired Legrendre polynomial. Next compute
            ** ppp its derivative by standard relation involving also p2,
            ** polynomial of one lower order.
            */
            pp = n * (z * p1 - p2)/(z * z - 1.0);
            z1 = z;
            z = z1 - p1/pp; // Newton's method
            //cout << "z = " << z << " z1 = " << z1 << endl;
        } while(fabs(z - z1) > ZERO);
        /*
        ** Scale the root to the desired interval and put in its symmetric
        ** counterpart. Compute the weight and its symmetric counterpart
        */
        *(x_low++) = xm - xl * z;
        *(x_high--) = xm + xl * z;
        *w_low = 2.0 * xl/((1.0 - z * z) * pp * pp);
        *(w_high--) = *(w_low++);
    }
} // End_ function gauleg()



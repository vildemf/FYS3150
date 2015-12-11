#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{
    //init. variables
    int n = 200;
    double delta_x = 5.0e-3;

    vec delta_t(4);
    delta_t(0) = 2.0e-8;
    delta_t(1) = 2.0e-7;
    delta_t(2) = 2.0e-6;
    delta_t(3) = 6.0e-6;

    for(int m = 0 ; m < 4 ; m++){

        double a = delta_t(m)/(delta_x*delta_x);

        //Make the matrix B
        mat I(n,n,fill::eye);
        mat B(n,n,fill::zeros);
        for(int i = 0 ; i < n ; i++){
            B(i,i) = 2;
            if(i<(n-1)){
                B(i,i+1) = -1;
                B(i+1,i) = -1;
            }
        }

        //Make the matrices for the three
        //methods; Bakcward Euler, Forward Euler
        //and Crank Nicolson scheme, and finding
        //their eigenvalues
        mat A_bkwrd(n,n);
        A_bkwrd = I + a*B;
        mat A_bkwrd_inv = inv(A_bkwrd);
        vec eigval_bkwrd = eig_sym(A_bkwrd_inv);
        //cout << eigval_bkrwd << endl;

        mat A_forw(n,n);
        A_forw = I - a*B;
        vec eigval_forw = eig_sym(A_forw);
        //cout << eigval_forw << endl;

        mat A_CrN(n,n);
        mat A_CrN1(n,n);
        mat A_CrN2(n,n);
        A_CrN1 = 2*I + a*B;
        A_CrN2 = 2*I - a*B;
        mat A_CrN1_inv = inv(A_CrN1);
        A_CrN = A_CrN1_inv*A_CrN2;
        vec eigval_CrN = eig_sym(A_CrN);
        //cout << eigval_CrN << endl;

        //Check if the absolute value of all
        //eigenvalues are bigger that one.
        for(int k = 0 ; k < n ; k++){
            double a = eigval_bkwrd(k);
            double b = eigval_forw(k);
            double c = eigval_CrN(k);
            if(abs(a) >= 1.0){
                cout << "Backward Euler "
                        "matrix not stable: eigval = "
                     << a << endl;
                exit(1);
            }
            if(abs(b) >= 1.0){
                cout << "Forward Euler "
                        "matrix not stable: eigval = "
                     << b << endl;
                exit(1);
            }
            if(abs(c) >= 1.0){
                cout << "Crank Nicolson "
                        "matrix not stable: eigval = "
                     << c << endl;
                exit(1);
            }

        }
    cout << "Tests pass for delta_t = " << delta_t(m) << endl;
    }

    return 0;
}



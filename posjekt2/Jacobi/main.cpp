#include <iostream>
#include <armadillo>
#include <cassert>
#include <fstream>

using namespace std;
using namespace arma;

// Declare functions before main()
double potential1(double rho);
double potential2(double rho, double omega);

// Jacobi functions
void jacobi_method(mat &A, mat &R, int n);
double find_max_offdiag(const mat &A, int &k, int &l, int n);
void rotate(mat &A, mat &R, int l, int k, int iterations, int n);
void compute_c_s(const mat &A, double &c, double &s, int k, int l);
void compute_A(mat &A, mat &R, double c, double s, int k, int l, int iterations, int n);
void compute_R(mat &R, double c, double s, int k, int l, int i);

// Test functions
void test_jacobi_method_2x2();
void test_compute_A();
void test_find_max_offdiag();
void test_compute_R_orthogonal(mat &R, int n); // Called at the end in jacobi_method()

int main()
{
    int N_step = 200;
    double rho_min = 0;
    double rho_max = 8;
    double h = (rho_max - rho_min)/N_step;

    // Set up rho, diagonal and off diagonal elements in vectors:
    // i=0 and i=N_step solutions known.
    // want indices 1-(N_step-1), that is 0-(N_step-2) in c++ indices.
    int n = N_step - 2;
    vec rho_list = linspace<vec>(rho_min + h, rho_max - h, n);
    double offdiag = -1.0/(h*h);
    // vec offdiag_list(n - 1); // not needed as vector..
    // offdiag_list.fill(offdiag);
    vec diag_list(n);
    double omega = 1.5;
    for (int i=0; i<n; i++) {
        diag_list(i) = 2.0/(h*h) + potential2(rho_list(i), omega);
    }
    mat A(n,n);
    // Physics ready. Set up matrix
    for (int i=0; i<n; i++) {
        A(i,i) = diag_list(i);
        if (i!=(n-1)) { // offdiagonals are 1 element fewer
            A(i,i+1) = offdiag;
            A(i+1,i) = offdiag;
        }
    }
    mat A_comparing = A; // Store the original to use in armadillo exact

    mat R(n,n);
    // Compute!
    //jacobi_method(A, R, n);
    //vec eigval_computed = sort(A.diag());

    // Exact from Armadillo
    vec eigval_exact;
    mat eigvec_exact;
    eig_sym(eigval_exact, eigvec_exact, A_comparing);
    //for (int i=0; i<n; i++) {
    //    double diff = abs(eigval_exact(i)-eigval_computed(i));
        //if (diff>1.0e-5) {
    //    cout << eigval_exact(i) << " " << eigval_computed(i) << endl;
        //}
    //}
    //cout << eigval_exact(0) << " " << eigval_computed(0) << endl;
    //cout << eigval_exact(1) << " " << eigval_computed(1) << endl;
    //cout << eigval_exact(2) << " " << eigval_computed(2) << endl;

    ofstream outfile0;
    outfile0.open ("/home/vilde/Documents/FYS3150/posjekt2/E0.txt");
    for (int i=0; i<n; i++) {
        outfile0 << eigvec_exact(i,0) << "\n";
    }
    outfile0.close();

    ofstream outfile1;
    outfile1.open ("/home/vilde/Documents/FYS3150/posjekt2/E1.txt");
    for (int i=0; i<n; i++) {
        outfile1 << eigvec_exact(i,1) << "\n";
    }
    outfile1.close();

    ofstream outfile2;
    outfile2.open ("/home/vilde/Documents/FYS3150/posjekt2/E2.txt");
    for (int i=0; i<n; i++) {
        outfile2 << eigvec_exact(i,2) << "\n";
    }
    outfile2.close();

    // Jacobi method testing
    test_compute_A();
    test_find_max_offdiag();
    test_jacobi_method_2x2();
    return 0;
}

double potential1(double rho) {
    return rho*rho;
}
double potential2(double rho, double omega) {
    return omega*omega*rho*rho + 1.0/rho;
}


void jacobi_method(mat &A, mat &R, int n) {
    int k = 0;
    int l = 0;
    R = eye<mat>(n,n);
    double epsilon = 1.0e-8;
    int max_iterations = n*n*n; // Set differently?
    int iterations = 0; // Iteration counter
    double max_offdiag = find_max_offdiag(A, k, l, n);

    // Square the max element in the test, as suggested in the text
    while (abs(max_offdiag) > epsilon && iterations < max_iterations) {
        rotate(A, R, l, k, iterations, n);
        max_offdiag = find_max_offdiag(A, k, l, n);
        iterations++;
    }
    if (iterations >= max_iterations) {
        cout << "OBS: Jacobi reached max iterations " << iterations << endl;
    }
    cout << "Iterations: " << iterations << endl;
    test_compute_R_orthogonal(R, n);

}


double find_max_offdiag(const mat &A, int &k, int &l, int n) {
    double max = 0.0;
    for (int i=0; i<n; i++) {
        for (int j=i+1; j<n; j++) { // Use that it is symmetrical
            if (abs(A(j,i)) > max) {
                max = abs(A(j,i));
                k = j;
                l = i;
            }
        }
    }
    return max;
}

void rotate(mat &A, mat &R, int l, int k, int iterations, int n) {
    double c = 0;
    double s = 0;
    compute_c_s(A, c, s, k, l);
    compute_A(A, R, c, s, k, l, iterations, n);
}

void compute_c_s(const mat &A, double &c, double &s, int k, int l) {
    if (A(k,l) == 0.0) {
        c = 1.0;
        s = 0.0;
    } else {
        double t = 0;
        double tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if (tau < 0) {
            t = -tau - sqrt(1 + tau*tau);
        } else {
            t = -tau + sqrt(1 + tau*tau);
        }
        c = 1.0/sqrt(1 + t*t);
        s = t*c;
    }
}

void compute_A(mat &A, mat &R, double c, double s, int k, int l, int iterations, int n) {
    // Store elements from old A needed to compute new A
    double a_kk = A(k,k); // reference, like A? How to keep from changing? const declaration, copy, or?
    double a_ll = A(l,l); // No it's fine they won't change (think this _is_ a copy. Not reference.).
    double a_ik = 0;
    double a_il = 0;
    double is_it_zero = (a_kk - a_ll)*c*s + A(k,l)*(c*c - s*s); // New a_kl should be zero by def. of c, s
    // Will it always be exactly zero? If not, should I set it to the small error number or hard code it as
    // zero, like I do now?
    //if (is_it_zero!=0.0) {
    //    cout << "OBS! The new a_kl is not zero at iteration " << iterations << " , it is" << endl; //
    //    cout << is_it_zero << endl;
    //}

    // Set new
    A(k,k) = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s;
    A(l,l) = a_ll*c*c + 2*A(k,l)*c*s + a_kk*s*s;
    A(k,l) = 0; // By def. of how we chose c,s, if everything went well (test above)
    A(l,k) = 0;

    // Set rest of new
    for (int i = 0; i<n; i++) {
        if (i!=l && i!=k) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = a_ik*c - a_il*s;
            A(k,i) = a_ik*c - a_il*s;
            A(i,l) = a_il*c + a_ik*s;
            A(l,i) = a_il*c + a_ik*s;
        }
        compute_R(R, c, s, k, l, i);
    }
    return;
}

void compute_R(mat &R, double c, double s, int k, int l, int i) {
    double r_ik = R(i,k);
    double r_il = R(i,l);
    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
}


void test_jacobi_method_2x2() {
    // Test matrix:
    int n = 2;
    mat B(n,n);
    B << 1 << 2 << endr << 2 << 4 << endr;
    mat R(n,n);
    // Solution from Matlab:
    mat B_exact(n,n);
    B_exact << 0.0 << 0.0 << endr << 0.0 << 5.0 << endr;
    jacobi_method(B, R, n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            // (j,i) not (i,j) because Armadillo says "mat"
            // is column major order.. not sure though
            assert(B(j,i) == B_exact(j,i));
        }
    }
}

void test_find_max_offdiag() {
    int n = 3;
    int k = 0;
    int l = 0;
    // Example real sym. matrix
    mat B(n,n);
    B << 4 << -2 << -50 << endr << -2 << -3 << 1 << endr << -50 << 1 << 1 << endr;
    double max = find_max_offdiag(B, k, l, n);
    assert(k==2 && l==0 && max==50);
}

void test_compute_A() {
    // Test matrix:
    int n = 2;
    int k = 0;
    int l = 1;
    mat R(n,n);
    mat B(n,n);
    B << 1 << 2 << endr << 2 << 4 << endr;
    double tau = 3.0/4;
    double t = -tau + sqrt(1 + tau*tau);
    double c = 1.0/sqrt(1 + t*t);
    double s = c*t;

    // Solution from Matlab:
    mat B_exact(n,n);
    B_exact << 0.0 << 0.0 << endr << 0.0 << 5.0 << endr;
    compute_A(B, R, c, s, k, l, 10, n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            // (j,i) not (i,j) because Armadillo says "mat"
            // is column major order.. not sure though
            assert(B(j,i) == B_exact(j,i));
        }
    }
}

void test_compute_R_orthogonal(mat &R, int n) {
    double tol = 1.0e-3;
    for (int i=0; i<n; i++) {
        vec eig1 = R.col(i);
        for (int j=0; j<n; j++) {
            //cout << R(i,j) << endl;
            vec eig2 = R.col(j);
            if (j!=i) { // i=j case would be 1. We test against 0. What if degenerate?
                double diff = abs(0 - dot(eig1,eig2));
                //cout << diff << endl;
                assert(diff < tol);
            }

        }
    }

}

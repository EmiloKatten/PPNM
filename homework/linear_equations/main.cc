#include <iostream>
#include <cmath>
#include <random>
#include "matrix.h"
#include "qr.h"

int main() {

	std::cout<< "======== TESTING QR CLASS========"<<"\n";
    // ----- 1. Generate random tall matrix A -----
    int n = 6;
    int m = 3;

    pp::matrix A(n, m);

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            A(i,j) = dist(gen);
    A.print("Random matrix A:");

    // ----- 2. QR decomposition -----
    pp::qr decomp(A);

    pp::matrix& Q = decomp.Q;
    pp::matrix& R = decomp.R;

    Q.print("Q:");
    R.print("R:");

    // ----- 3. Check R is upper triangular -----
    bool upper = true;

    for (int i = 0; i < R.size1(); ++i)
        for (int j = 0; j < i; ++j)
            if (std::abs(R(i,j)) > 1e-10)
                upper = false;

    std::cout << "R upper triangular: "
              << (upper ? "TRUE\n" : "FALSE\n");

    // ----- 4. Check Q^T Q = I -----
    pp::matrix Qt = Q.transpose();
    pp::matrix QtQ = Qt * Q;
	QtQ.threshold();
	QtQ.print("Q^T Q = I ?");

    // ----- 5. Check QR = A -----
    pp::matrix QR = Q * R;
    bool check_QR = pp::approx(QR, A);
    std::cout<<"Is QR = A? " << (check_QR ? "TRUE\n" : "FALSE\n");



	std::cout<<"\n\n\n"<<"======== TESTING QR SOLVE ========"<<"\n";
	
	// ----- 1. Generate random square matrix -----
	pp::matrix B(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            B(i,j) = dist(gen);
    B.print("Random square matrix B =");

	// ----- 2. Random vector -----
	pp::vector b(n);
    for (int i = 0; i < n; i++){
        b[i] = dist(gen);
	}
	b.print("Random vector b =");

	// ----- 3. Factorize B -----
	pp::qr decomp2(B);

	// ----- 4. Solve QRx = b -----
	pp::vector x = decomp2.solve(b);
	x.print("Solving QRx = b: \nx = ");

	// ----- 5. Check Bx = b -----
	pp::vector Bx(n);
	for (int i=0; i < n; i++){
		Bx[i] = pp::dot(B.transpose()[i], x);
	}
	std::cout<<"Is Bx = b? ";
	pp::vector check = Bx - b;
    bool test = true;
    for (int i = 0; i < n; i++){
        if (check[i] > 1e-10) {test = false;}
    }
    std::cout<< (test ? "TRUE" : "FALSE") <<"\n";


	std::cout<<"\n\n\n"<<"======== DETERMINANT OF R ========"<<"\n";
    R.print("R:");
    std::cout<<"Det. of R: "<<decomp.det()<<"\n";



	std::cout<<"\n\n\n"<<"======== INVERSE OF MATRIX ========"<<"\n";
    B.print("Random square matrix B:");
    pp::matrix B_inv = decomp2.inverse();
    B_inv.print("Inverse of matrix B:");

    pp::matrix BB_inv = B * B_inv;
    BB_inv.threshold();
    BB_inv.print("B * B_inv = I?");

    return 0;
}
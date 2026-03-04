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
	pp::matrix diff = QR - A;
	diff.threshold();
	diff.print("Is QR = A, then QR - A = 0:");



	std::cout<<"\n\n\n"<<"======== TESTING QR SOLVE ========"<<"\n";
	
	// ----- 1. Generate random square matrix -----
	pp::matrix A2(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A2(i,j) = dist(gen);
    A2.print("Random square matrix B =");

	// ----- 2. Random vector -----
	pp::vector b(n);
    for (int i = 0; i < n; i++){
        b[i] = dist(gen);
	}
	b.print("Random vector b =");

	// ----- 3. Factorize A -----
	pp::qr decomp2(A2);

	// ----- 4. Solve QRx = b -----
	pp::vector x = decomp2.solve(b);
	x.print("Solving QRx = b: \nx = ");

	// ----- 5. Check Ax = b -----
	pp::vector Ax(n);
	for (int i=0; i < n; i++){
		Ax[i] = pp::dot(A2.transpose()[i], x);
	}
	std::cout<<"Is Bx = b? ";
	pp::vector check = Ax - b;
    bool test = true;
    for (int i = 0; i < n; i++){
        if (check[i] > 1e-10) {test = false;}
    }
    std::cout<< (test ? "TRUE" : "FALSE") <<"\n";


	std::cout<<"\n\n\n"<<"======== DETERMINANT OF R ========"<<"\n";
    R.print("R:");
    std::cout<<"Det. of R: "<<decomp.det()<<"\n";



	std::cout<<"\n\n\n"<<"======== INVERSE OF MATRIX ========"<<"\n";
    A2.print("Random square matrix B:");
    pp::matrix A_inv = decomp2.inverse();
    A_inv.print("Inverse of matrix B:");

    pp::matrix AA_inv = A2 * A_inv;
    AA_inv.threshold();
    AA_inv.print("B * B_inv = I?");

    return 0;
}
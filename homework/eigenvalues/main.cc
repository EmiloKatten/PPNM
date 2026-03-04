#include<random>
#include<iostream>
#include"matrix.h"
#include"evd.h"

int main(){

    //----- Make random symmetric matrix A -----
    int n = 6;
    pp::matrix A(n, n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double value = dist(gen);
            A(i, j) = value;
            A(j, i) = value;  // mirror
        }
    }

    //----- Print A, V, w and D -----
    pp::EVD evdA(A);
    pp::matrix& V = evdA.V;
    pp::vector& w = evdA.w;
    A.print("Random symmetric matrix A");
    V.print("V:");
    w.print("w:");

    pp::matrix D(n, n);
    for (int i = 0; i < n; i++){
        D(i,i) = w[i];
    }
    D.print("D:");

    
    //----- Checks -----
    std::cout<<"\n\n\n=====TESTING====="<<"\n";

    pp::matrix A_test = V*D*V.transpose();
    bool A_check = true;
    for (int i = 0; i < n; i++){
        if (!pp::approx(A[i], A_test[i])) {A_check = false;}
    }
    A_test.print("VDV^T:");
    std::cout<<"Is A == VDV^T? "<<(A_check ? "TRUE\n" : "FALSE\n");


    pp::matrix D_test = V.transpose()*A*V;
    D_test.threshold();
    bool D_check = true;
    for (int i = 0; i < n; i++){
        if (!pp::approx(D[i], D_test[i])) {D_check = false;}
    }
    D_test.print("\nV^TAV:");
    std::cout<<"Is D == V^TAV? "<<(D_check ? "TRUE\n" : "FALSE\n");


    pp::matrix V_test1 = V*V.transpose();
    pp::matrix V_test2 = V.transpose()*V;
    V_test1.threshold();
    V_test2.threshold();

    V_test1.print("\nVV^T:");
    bool V_check = true;
    for (int i = 0; i < n; i++){
        if (!pp::approx(V_test1[i], V_test2[i])) {V_check = false;}
    }
    std::cout<<"Is VV^T == V^TV? "<<(V_check ? "TRUE\n" : "FALSE\n");


    return 0;
}
#include<random>
#include<iostream>
#include<fstream>
#include"matrix.h"
#include"evd.h"

pp::vector make_r(double rmax, double dr){
    int npoints = (int)(rmax/dr)-1;
    pp::vector r(npoints);
    for(int i=0;i<npoints;i++){r[i]=dr*(i+1);}
    return r;
}

pp::matrix make_H(double rmax, double dr){
    pp::vector r = make_r(rmax, dr);
    int npoints=r.size();
    pp::matrix H(npoints,npoints);
    for(int i=0;i<npoints-1;i++){
        H[i,i]  =-2*(-0.5/dr/dr);
        H[i,i+1]= 1*(-0.5/dr/dr);
        H[i+1,i]= 1*(-0.5/dr/dr);
    }
    H[npoints-1,npoints-1]=-2*(-0.5/dr/dr);
    for(int i=0;i<npoints;i++)H[i,i]+=-1/r[i];
    return H;
}

double compute_epsilon0(double rmax, double dr){
    pp::matrix H = make_H(rmax, dr);
    pp::EVD evdH(H);
    pp::vector& w = evdH.w;
    double epsilon0 = w[0];
    for (int i = 1; i < w.size(); i++){
        if (w[i] < epsilon0) {epsilon0 = w[i];}
    }
    return epsilon0;
}


int main(int argc,char** argv){

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



    //=======EXERCISE B=========
    double rmax = 0.0;
    double dr   = 0.0;

    // Loop through command-line arguments
    for(int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if(arg == "-rmax" && i+1 < argc) {
            rmax = std::stod(argv[++i]);
        }
        else if(arg == "-dr" && i+1 < argc) {
            dr = std::stod(argv[++i]);
        }
    }

    std::cout<<"\n\n\n===== BUILDING HAMILTONIAN MATRIX =====\n";
    if (rmax == 0 || dr == 0){
        std::cout<<"No arguments given - skipping\n";
        return 0;
    }
    std::cout<<"-rmax = "<<rmax<<" -dr = "<<dr<<"\n";
    double eps0 = compute_epsilon0(rmax, dr);
    std::cout << "Epsilon_0 (ground energy): " << eps0 << "\n";

    // Construct data for fixed rmax and dr
    std::ofstream data_dr("data_dr.txt");
    for (double i = 0.05; i < 1; i+=0.05){
        double eps0_dr = compute_epsilon0(rmax, i);
        data_dr<<i<<" "<<eps0_dr<<"\n";
    }
    std::ofstream data_rmax("data_rmax.txt");
    for (double i = 1; i < 10; i+=0.3){
        double eps0_rmax = compute_epsilon0(i, dr);
        data_rmax<<i<<" "<<eps0_rmax<<"\n";
    }


    // Construct lowest order eigen-functions and write data
    pp::matrix H = make_H(rmax, dr);
    pp::vector r = make_r(rmax, dr);

    int npoints=r.size();
    int k_values=3;

    pp::EVD evdH(H);
    pp::matrix& V_H = evdH.V;
    pp::vector& w_H = evdH.w;

    // index array
    int* idx=new int[npoints];
    for(int i=0;i<npoints;i++) idx[i]=i;

    // selection sort by eigenvalue
    for(int i=0;i<npoints-1;i++){
        int min=i;
        for(int j=i+1;j<npoints;j++)
            if(w_H[idx[j]]<w_H[idx[min]]) min=j;
        std::swap(idx[i],idx[min]);
    }

    std::ofstream file("eigenfunctions.txt");
    file<<"0 "<<"0 "<<"0 "<<"0 "<<"0 "<<"0 "<<"0\n";
    double Const=1/std::sqrt(dr);
    for(int i=0;i<npoints;i++){
        file<<r[i];
        for(int k=0;k<k_values;k++){
            double f=Const*V_H(i,idx[k]);
            file<<" "<<f;
        }
        double factor = V_H[0][0]>0 ? 1:-1;
        double f1 = 2*r[i]*std::exp(-r[i])*factor;
        file<<" "<<f1;
        factor = V_H[1][0]>0 ? 1:-1;
        double f2 = 1/std::sqrt(2)*r[i]*(1-r[i]/2.0)
                    * std::exp(-r[i]/2.0)*factor;
        file<<" "<<f2;
        factor = V_H[2][0]>0 ? 1:-1;
        double f3 = 2.0/3.0/std::sqrt(3)*r[i]
                    * (1-2*r[i]/3.0+2*r[i]*r[i]/27.0)
                    * std::exp(-r[i]/3.0)*factor;
        file<<" "<<f3;

        file<<"\n";
    }



    return 0;
}
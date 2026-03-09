#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <functional>
#include <fstream>
#include "matrix.h"
#include "qr.h"
#include "ls.h"


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
    A.print("Random tall matrix A:");

    // ----- 2. QR decomposition -----
    pp::qr decomp(A);

    pp::matrix& Q = decomp.Q;
    pp::matrix& R = decomp.R;

    Q.print("Q:");
    R.print("R:");


    std::cout<<"\n\n\n======= FITTING TO DATA USING QR DECOMPOSITION =======\n";

    std::vector<std::function<double(double)>> fs {
        [](double z) { return 1.0; }, //ln(a)
        [](double z) { return z; }, //-lambda
    };

    // Writing data
    pp::vector t = {1, 2, 3, 4, 6, 9, 10, 13, 15};
    pp::vector y = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};
    pp::vector lny(y.size());
    pp::vector dy = {6, 5, 4, 4, 4, 3, 3, 2, 2};
    pp::vector dlny(dy.size());
    for (int i=0; i < y.size(); i++){
        lny[i] = std::log(y[i]);    // taking ln(y)
        dlny[i] = dy[i] / y[i];     // computing delta ln(y)
    }
    pp::vector fit = pp::lsfit(fs, t, lny, dlny);   // getting fit parameters
    double lna = fit[0];
    double lambda = fit[1];
    std::cout<<"Fitting data to ln(y)=ln(a)-λt\n\tln(a) = "<<lna<<"\n\tλ = "<<-1*lambda<<"\n";
    pp::vector fit_lny(y.size());
    pp::vector fit_y(y.size());
    for (int i=0; i < y.size(); i++){
        fit_lny[i] = lna + lambda*t[i];
        fit_y[i] = std::exp(fit_lny[i]);
    }

    //  writing to data file
    std::ofstream file("data.txt");
    for (int i=0; i < t.size(); i++){
        file<<t[i]<<" "<<y[i]<<" "<<dy[i]<<" "<<fit_y[i]<<"\n";
    }
    

    // Compute half-life
    double half_life = -std::log(2)/lambda;
    std::cout<<"Half-life of fit: "<<half_life<<" days \n";
    // 224 Ra half-life: 3.6316(23) days
    std::cout<<"Half-life of 224Ra: 3.6316(23) days \n";
    return 0;
}
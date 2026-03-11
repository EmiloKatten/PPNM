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
    auto [fit, sigma] = pp::lsfit(fs, t, lny, dlny);   // getting fit parameters + uncertainties
    double lna = fit[0];
    double lambda = fit[1];
    std::cout<<"Fitting data to ln(y)=ln(a)-λt\n\tln(a) = "<<lna<<"\n\tλ = "<<-1*lambda<<"\n";

    //  writing to data file
    std::ofstream file_data("data.txt");
    for (int i=0; i < t.size(); i++){
        file_data<<t[i]<<" "<<y[i]<<" "<<dy[i]<<" "<<"\n";
    }
    std::ofstream file_fit("fit.txt");
    for (double i=t[0]; i < t[t.size()-1]; i+=0.1){ //smooth fit
        double fit_lny = lna + lambda*i;
        double fit_y = std::exp(fit_lny);
        file_fit<<i<<" "<<fit_y<<"\n";
    }
    

    // Compute half-life
    double half_life = -1*std::log(2)/lambda;
    std::cout<<"Half-life of fit: "<<half_life<<" days \n";
    // 224 Ra half-life: 3.6316(23) days
    std::cout<<"Half-life of 224Ra: 3.6316(23) days \n";



    std::cout<<"\n\n\n======= UNCERTAINTIES OF THE FITTING PARAMETERS =======\n";
    sigma.print("Uncertainty matrix:");
    double lna_err = std::sqrt(sigma[0][0]);
    double lambda_err = std::sqrt(sigma[1][1]);
    std::cout<<"\t ln(a) error = "<<lna_err<<"\n";
    std::cout<<"\t λ error = "<<lambda_err<<"\n";

    //Using error propagation, T_1/2_err = |d/dλ T_1/2| * λ_err = ln(2)/λ^2 * λ_err, so:
    double half_life_err = std::log(2)/lambda/lambda * lambda_err;
    std::cout<<"\nUsing error propagation:\n\tHalf-life error = "<<half_life_err<<"\n";
    std::cout<<"\nThus, fit half-life: "<<half_life<<" ± "<<half_life_err<<"\n";
    std::cout<<"It does not quite agree with the modern value... ☹";


    return 0;
}
#include "mc.h"
#include "lcg.h"
#include<iostream>
#include<fstream>

int main(){

    // Define function f(x,y)
    auto f = [](const vec& x){
        double r2 = x[0]*x[0] + x[1]*x[1];  //r^2 = x^2 + y^2
        return (r2 <= 1.0) ? 1.0 : 0.0;     //1 inside circle, 0 outside circle
    };

    // Integration limits
    vec a = {-1.0, -1.0};
    vec b = { 1.0,  1.0};
    pp::lcg rng(15042026);  // random number generator with seed

    auto [result, error] = pp::plainmc(f, a, b, 1e6, rng);

    std::cout << "EXERCISE A\n\nUnit circle integral...\n";
    std::cout << "\tEstimated area = " << result << "\n";
    std::cout << "\tEstimated error = " << error << "\n";
    std::cout << "\tActual value = " << M_PI << "\n";
    std::cout << "\tActual error = " << std::abs(result - M_PI) << "\n";



    // Plot errors:
    std::ofstream data_error("data_error.txt");
    for(int N = 10; N<1e9; N*=10){
        auto [results, errors] = pp::plainmc(f, a, b, N, rng);
        data_error << N << " " << errors << " " << std::abs(results - M_PI) << " " << 1/std::sqrt(N) << "\n";
    }


    // Compute difficult integral:
    auto g = [](const vec& x){
        double res = 1/std::pow(M_PI,3) / (1-std::cos(x[0])*std::cos(x[1])*std::cos(x[2]));
        return res;
    };
    vec a_g = {0.0, 0.0, 0.0};
    vec b_g = {M_PI, M_PI, M_PI};
    auto [result_g, error_g] = pp::plainmc(g, a_g, b_g, 1e6, rng);
    double actual_result_g = 1.3932039296856768591842462603255;

    std::cout << "\nDifficult integral...\n"; 
    std::cout << "\tEstimated value = " << result_g << "\n";
    std::cout << "\tEstimated error = " << error_g << "\n";
    std::cout << "\tActual value = " << actual_result_g << "\n";
    std::cout << "\tActual error = " << std::abs(result_g - actual_result_g) << "\n";



    // ==== PART B ====
    auto qresults = pp::quasimc(f, a, b, 1e6);
    std::cout << "\n\nEXERCISE B\n\nUnit circle integral...\n";
    std::cout << "\tEstimated area = " << qresults << "\n";
    std::cout << "\tActual value = " << M_PI << "\n";
    std::cout << "\tActual error = " << std::abs(qresults - M_PI) << "\n";

    std::ofstream quasi_error("quasi_error.txt");
    for(int N = 10; N<1e9; N*=10){
        auto qresults_ = pp::quasimc(f, a, b, N);
        quasi_error << N << " " << std::abs(qresults_ - M_PI) << "\n";
    }

    return 0;
}
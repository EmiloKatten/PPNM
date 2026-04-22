#include "mc.h"
#include "lcg.h"
#include<iostream>
#include<fstream>
#include <random>

int main(){

    // Define function f(x,y)
    auto f_circle = [](const vec& x){
        double r2 = x[0]*x[0] + x[1]*x[1];  //r^2 = x^2 + y^2
        return (r2 <= 1.0) ? 1.0 : 0.0;     //1 inside circle, 0 outside circle
    };

    // Integration limits
    vec a_circle = {-1.0, -1.0};
    vec b_circle = { 1.0,  1.0};
    pp::lcg rng(15042026);  // random number generator with seed

    auto [result_circle, error_circle] = pp::plainmc(f_circle, a_circle, b_circle, 1e6, rng);

    std::cout << "EXERCISE A\n\nUnit circle integral...\n";
    std::cout << "\tEstimated area = " << result_circle << "\n";
    std::cout << "\tEstimated error = " << error_circle << "\n";
    std::cout << "\tActual value = " << M_PI << "\n";
    std::cout << "\tActual error = " << std::abs(result_circle - M_PI) << "\n";


    // 3D ellipsoid
    auto f_ellipsoid = [](const vec& x){
        double val =
            (x[0]*x[0])/(1.0*1.0) +   // x^2 / a^2
            (x[1]*x[1])/(2.0*2.0) +   // y^2 / b^2
            (x[2]*x[2])/(3.0*3.0);    // z^2 / c^2
        return (val <= 1.0) ? 1.0 : 0.0;
    };
    vec a_ellipsoid = {-1.0, -2.0, -3.0};
    vec b_ellipsoid = { 1.0,  2.0,  3.0};
    auto [result_ellipsoid, error_ellipsoid] = pp::plainmc(f_ellipsoid, a_ellipsoid, b_ellipsoid, 1e6, rng);
    double true_volume = (4.0/3.0)*M_PI*1*2*3;

    std::cout << "3D ellipsoid...\n";
    std::cout << "\tEstimated volume = " << result_ellipsoid << "\n";
    std::cout << "\tEstimated error = " << error_ellipsoid << "\n";
    std::cout << "\tActual value = " << true_volume << "\n";
    std::cout << "\tActual error = " << std::abs(result_ellipsoid - true_volume) << "\n";

    // Plot errors:
    std::ofstream data_error("data_error.txt");
    for(int N = 10; N<1e9; N*=10){
        auto [results, errors] = pp::plainmc(f_circle, a_circle, b_circle, N, rng);
        data_error << N << " " << errors << " " << std::abs(results - M_PI) << " " << 1/std::sqrt(N) << "\n";
    }


    
    // ==== PART B ====
    // Compute difficult integral:
    auto g = [](const vec& x){
        double res = 1/std::pow(M_PI,3) / (1-std::cos(x[0])*std::cos(x[1])*std::cos(x[2]));
        return res;
    };
    vec a_g = {0.0, 0.0, 0.0};
    vec b_g = {M_PI, M_PI, M_PI};
    auto [result_g_lcg, error_g_lcg] = pp::plainmc(g, a_g, b_g, 1e6, rng);
    std::mt19937 rng_std(15042026);  // seed
    auto [result_g_std, error_g_std] = pp::plainmc(g, a_g, b_g, 1e6, rng_std);
    auto qresults = pp::quasimc(g, a_g, b_g, 1e6);  
    double actual_result_g = 1.3932039296856768591842462603255;

    std::cout << "\n\nEXERCISE B\n\nDifficult integral...\n"; 
    std::cout << "\tActual value = " << actual_result_g << "\n";
    std::cout << "\tEstimated value (LCG)= " << result_g_lcg << " error = "<< std::abs(result_g_lcg - actual_result_g) <<"\n";
    std::cout << "\tEstimated value (STD)= " << result_g_std << " error = "<< std::abs(result_g_std - actual_result_g) <<"\n";
    std::cout << "\tEstimated value (quasi)= " << qresults << ", error = "<< std::abs(qresults - actual_result_g) <<"\n";

    return 0;
}
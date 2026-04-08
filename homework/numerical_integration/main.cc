#include "integrator.h"
#include <iostream>
#include <cmath>
#include <fstream>

int main(){
    
    // start by testing the integration:
    int ncalls1 = 0;
    auto f1 = [&ncalls1](double x){
        ncalls1++;
        return std::sqrt(x);
    };
    int ncalls2 = 0;
    auto f2 = [&ncalls2](double x){
        ncalls2++;
        return 1.0 / std::sqrt(x);
    };
    int ncalls3 = 0;
    auto f3 = [&ncalls3](double x){
        ncalls3++;
        return std::sqrt(1-x*x);
    };
    int ncalls4 = 0;    
    auto f4 = [&ncalls4](double x){
        ncalls4++;
        return std::log(x) / std::sqrt(x);
    };
    double q1 = pp::integrate(f1,0,1); bool approx1 = pp::approx(q1, 2.0/3.0);
    double q2 = pp::integrate(f2,0,1); bool approx2 = pp::approx(q2, 2.0);
    double q3 = pp::integrate(f3,0,1); bool approx3 = pp::approx(q3, pp::PI/4.0);
    double q4 = pp::integrate(f4,0,1); bool approx4 = pp::approx(q4, -4.0);
    std::cout << std::boolalpha;
    std::cout << "sqrt(x): q="<<q1<<".   Within accuracy: "<<approx1<<".   ncalls="<<ncalls1<<"\n";
    std::cout << "1/sqrt(x): q="<<q2<<".   Within accuracy: "<<approx2<<".   ncalls="<<ncalls2<<"\n";
    std::cout << "sqrt(1-x^2): q="<<q3<<".   Within accuracy: "<<approx3<<".   ncalls="<<ncalls3<<"\n";
    std::cout << "log(x)/sqrt(x): q="<<q4<<".   Within accuracy: "<<approx4<<".   ncalls="<<ncalls4<<"\n";


    // ### TABULATING ERF ###
    std::ofstream data_erf("data_erf.txt");
    for(double i = -2.0; i <= 2.0; i+=0.1){
        data_erf<< i <<" "<< pp::erf(i) <<"\n";
    } 
    double const erf1 = 0.84270079294971486934;
    std::ofstream data_erf_acc("data_erf_acc.txt");
    for(double acc = 0.1; acc >= 1e-9; acc/=10){
        double tab_erf1 = pp::erf(1, acc, 0.0);
        data_erf_acc<< acc <<" "<< std::abs(erf1-tab_erf1) <<"\n";
    }



    // ######## EXERCISE B #########
    std::cout << "\nNow using variable transformation (Clenshaw-Curtis)\n";
    int cc_ncalls1 = 0;
    auto cc_f1 = [&cc_ncalls1](double x){
        cc_ncalls1++;
        return 1.0 / std::sqrt(x);
    };
    int cc_ncalls2 = 0;
    auto cc_f2 = [&cc_ncalls2](double x){
        cc_ncalls2++;
        return std::log(x) / std::sqrt(x);
    };
    double cc_q1 = pp::cc_integrate(cc_f1, 0, 1);
    double cc_q2 = pp::cc_integrate(cc_f2, 0, 1);
    std::cout << "1/sqrt(x): q="<<cc_q1<<".   ordinary ncalls="<<ncalls2<<",   variable transformation ncalls="<<cc_ncalls1<<",   python ncalls=231"<<"\n";
    std::cout << "log(x)/sqrt(x): q="<<cc_q2<<".   ordinary ncalls="<<ncalls4<<",   variable transformation ncalls="<<cc_ncalls2<<",   python ncalls=315"<<"\n";
    
    
    // ### TEST (CONVERGING) INFINITE INTEGRAL ###
    int cc_ncalls3 = 0;
    auto cc_f3 = [&cc_ncalls3](double x){
        cc_ncalls3++;
        return std::exp(-x);
    };
    double cc_q3 = pp::integrate(cc_f3, 0, pp::inf, 1e-4, 1e-4); bool cc_approx = pp::approx(cc_q3, 1.0);
    std::cout << "\n Testing infinite integral.\n";
    std::cout << "exp(-x) from 0 to inf: q="<<cc_q3<<".   Within accuracy: "<<cc_approx<<".   ncalls="<<cc_ncalls3<<",   python ncalls=135\n";

    return 0;
}

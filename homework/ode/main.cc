#include<tuple>
#include<functional>
#include<cmath>
#include<algorithm>
#include<iostream>
#include<numbers>
#include<fstream>
#include"matrix.h"
#include"rkstep.h"

double pi = std::numbers::pi;

using vector = pp::vector;


int main(){

    // u'' = -u
    // something like u = cos(x)
    // let y[0] = u and y[1] = u'
    std::function<vector(double,vector)> f1 = [](double x, vector y){
        vector dydx(2);
        dydx[0] = y[1];   // u' = dy[0]/dx = y[1]
        dydx[1] = -y[0];  // u'' = dy[1]/dx = -y[0]
        return dydx;
    };

    double start1 = 0;   //start point
    double end1 = pi;  //end point

    vector cond1(2);
    //start conditions
    cond1[0] = 1;  // cos(0) = 1
    cond1[1] = 0;  // cos'(0) = -sin(0) = 0

    auto [xlist1, ylist1] = pp::driver(f1, start1, end1, cond1);

    // Print results
    std::cout<<"Computing ODE: u'' = -u \n";
    std::cout<<"cos(π) = "<<std::cos(end1)<<", rkstep12 at π: "<<ylist1[xlist1.size()-1][0]<<"\n";
 


    // ========= OSCILLATOR WITH FRICTION =========
    double b = 0.25; // damping coefficient
    double c = 5.0; // gravity / length factor

    std::function<vector(double, vector)> f2 = [b,c](double t, vector y) {
        vector dydt(2);
        dydt[0] = y[1];                       // θ' = y1
        dydt[1] = -b*y[1] - c*std::sin(y[0]); // θ'' = -b*θ' - c*sin(θ)
        return dydt;
    };

    double start2 = 0;
    double end2 = 10;

    vector cond2(2);
    cond2[0] = pi - 0.1;
    cond2[1] = 0.0;

    auto [xlist2, ylist2] = pp::driver(f2, start2, end2, cond2);

    // Save results
    std::ofstream data_oscillator("data_oscillator.txt");
    for(size_t i = 0; i < xlist2.size(); i++){
        data_oscillator<< xlist2[i] <<" "<< ylist2[i][0] <<" "<< ylist2[i][1]<<"\n";
    }



    // ======= EXERCISE B =======
    double eps = 0.01;  // damping coefficient

    std::function<vector(double, vector)> f3 = [](double x, vector y) {
        vector dydx(2);
        dydx[0] = y[1];
        dydx[1] = 1-y[0];
        return dydx;
    };
    std::function<vector(double, vector)> f3eps = [eps](double x, vector y) {
        vector dydx(2);
        dydx[0] = y[1];
        dydx[1] = 1-y[0]+eps*y[0]*y[0];
        return dydx;
    };

    double start3 = 0;
    double end3 = 10*pi;

    vector cond3(2);
    cond3[0] = 1.0;
    cond3[1] = 0.0;

    vector cond3eps(2);
    cond3eps[0] = 1.0;
    cond3eps[1] = -0.5;

    double h=0.1,acc=0.01,relerr=0.01,hmin=5.0/360;
    auto [xlist3, ylist3] = pp::driver(f3, start3, end3, cond3, h, acc, relerr, hmin);
    auto [xlist4, ylist4] = pp::driver(f3, start3, end3, cond3eps);
    auto [xlist5, ylist5] = pp::driver(f3eps, start3, end3, cond3eps);


    // Save results
    std::ofstream data_planets1("data_planets1.txt");
    for(size_t i = 0; i < xlist3.size(); i++){
        data_planets1<< xlist3[i] <<" "<< ylist3[i][0] <<"\n";
    }
    std::ofstream data_planets2("data_planets2.txt");
    for(size_t i = 0; i < xlist4.size(); i++){
        data_planets2<< xlist4[i] <<" "<< ylist4[i][0] <<"\n";
    }
    std::ofstream data_planets3("data_planets3.txt");
    for(size_t i = 0; i < xlist5.size(); i++){
        data_planets3<< xlist5[i] <<" "<< ylist5[i][0] <<"\n";
    }


    return 0;
}

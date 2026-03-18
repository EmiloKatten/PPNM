#include<vector>
#include<cassert>
#include<cmath>
#include<fstream>
#include<iostream>
#include"qspline.h"

using vec = std::vector<double>;

int main(){

    // MAKING LINEAR INTERPOLATION
    std::cout<<"making linterp"<<"\n";
    vec data_x; vec data_y;
    std::ofstream file_data("ldata.txt");
    for (int i=0; i < 19; i++){
        data_x.push_back(0.5*i);
        data_y.push_back(std::cos(data_x[i]));
        file_data<<data_x[i]<<" "<<data_y[i]<<"\n";
    }
    std::cout<<"writing linterp"<<"\n";
    std::ofstream file_linterp("linterp.txt");
    for (double z=0; z<9; z+=0.1){
        double interp = pp::linterp(data_x, data_y, z);
        double integ = pp::linterpInteg(data_x, data_y, z);
        file_linterp<<z<<" "<<interp<<" "<<integ<<"\n";
    }


    // MAKING QUADRATIC INTERPOLATION
    double n_points = 10.0;
    std::cout<<"making qinterp"<<"\n";
    vec quad_x, quad_y, quad_y2;
    std::ofstream file_qdata("qdata.txt");
    for (double x=1; x<5.5; x+=5/n_points){
        double y = x*x;
        double y2 = std::cos(3*x)*x*x+10;

        quad_x.push_back(x);
        quad_y.push_back(y);
        quad_y2.push_back(y2);

        file_qdata<<x<<" "<<y<<" "<<y2<<"\n";
    }
    std::cout<<"writing qinterp"<<"\n";
    pp::qspline splinef(quad_x, quad_y);
    pp::qspline splineg(quad_x, quad_y2);
    std::ofstream file_qinterp("qinterp.txt");
    for (double z = 1; z<5; z+=0.1){
        double f = splinef.eval(z);
        double df = splinef.deriv(z);
        double F = splinef.integ(z);
        double g = splineg.eval(z);

        double ff = z*z;
        double dff = 2*z;
        double FF = 1/3.0 * z*z*z;
        double gg = std::cos(3*z)*z*z+10;
        file_qinterp<<z<<" "<<f<<" "<<ff<<" "<<df<<" "<<dff<<" "<<F<<" "<<FF<<" "<<g<<" "<<gg<<"\n";
    }
}
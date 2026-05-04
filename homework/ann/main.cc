#include<iostream>
#include<fstream>
#include<cmath>
#include"ann.h"

int main(){

    std::cout<<"EXERCISE A\n\n";
    // sample som xs and ys for ann //
    std::ofstream out("data.txt");
    pp::vector xs; pp::vector ys; 
    for (double x=-1; x<1; x+=0.05){
        double fx = std::cos(5*x-1)*std::exp(-x*x);
        xs.push_back(x); ys.push_back(fx);
        out << x << " " << fx << "\n";
    }

    pp::ann func(5);
    func.train(xs, ys);
    //testing:
    std::cout<<"g(x)=cos(5x-1)*exp(-x^2)\nFinished training\n";
    std::cout<<"\tg(0)="<<std::cos(-1.0)<<", f(0)="<<func.response(0, func.p)<<"\n";

    std::ofstream out2("fit.txt");
    for (double x = -1; x < 1; x += 0.05){
        double pred_y = func.response(x, func.p);
        out2 << x << " " << pred_y << "\n";
    }


    std::cout<<"\n\nEXERCISE B\n\n";
    // check derivatives and anti-derivative
    std::cout<<"Computing derivatives and anti-derivative\n";
    std::cout<<"\tg'(0.5)=-3.93934"<<", f'(0.5)="<<func.diff1(0.5, func.p)<<"\n";
    std::cout<<"\tg''(0.5)=6.33615"<<", f''(0.5)="<<func.diff2(0.5, func.p)<<"\n";
    std::cout<<"\tG(0.5)=0.158085"<<", F(0.5)="<<func.integral(0.5, func.p)<<"\n";
    // debugging
    std::ofstream out3("integral.txt");
    for (double x = -1; x < 1; x += 0.05){
        double pred_y = func.integral(x, func.p);
        out3 << x << " " << pred_y << "\n";
    }
    std::cout<<"The anti-derivative looks correct, but it's missing an offset factor\n";
    std::cout<<"This could have been implemented with a bias parameter\n";
    return 0;
}
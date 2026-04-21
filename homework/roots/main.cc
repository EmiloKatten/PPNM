#include"root.h"
#include"rkstep.h"
#include<iostream>
#include <fstream>

/* M function made with help from chatbot */
std::tuple<std::vector<double>, std::vector<pp::vector>> M(double E, double rmin, double rmax, double acc=0.01, double eps=0.01){
    std::function<pp::vector(double, pp::vector)> F = [E](double r, pp::vector y)
    {
        pp::vector dydr(2);

        double f  = y[0];
        double fp = y[1];

        dydr[0] = fp;
        dydr[1] = -2.0 * (1.0/r + E) * f;

        return dydr;
    };

    pp::vector y0 = {rmin - rmin*rmin, 1 - 2*rmin};

    return pp::driver(F, rmin, rmax, y0, 0.125, acc, eps);
};

double compute_E(pp::vector E0, double rmin = 1e-3, double rmax = 8, double acc = 0.01, double eps=0.01){ 
    std::function<pp::vector(pp::vector)> M_func = [&](pp::vector E0){
        double E = E0[0];
        pp::vector y(1);
        auto [rlist, ylist] = M(E, rmin, rmax, acc, eps); 
        y[0] = ylist.back()[0];
        return y;
    };
    pp::vector root_M = pp::newton(M_func, E0);
    return root_M[0];
};
    
double true_wavefunction(double r){
    return r * std::exp(-r);
};

int main(){

    std::cout<<"EXERCISE A\n\n";

    std::function<pp::vector(pp::vector)> debugf = [](pp::vector x) {
        pp::vector y(1);
        y[0] = x[0]*x[0] - 4;
        return y;
    };
    pp::vector x0(1); x0[0] = 1.0;  // inital guess, should be +-2
    pp::vector root_debugf = pp::newton(debugf, x0);
    std::cout<<"quick debugging: root of x^2-4: "<< root_debugf[0] <<"\n";


    //Rosenbrock's valley function
    // f(x,y) = (1-x)^2 + 100(y-x^2)^2
    // df/dx = 2*(200x^3 - 200xy + x-1)
    // df/dy = 200*(y - x^2)

    std::function<pp::vector(pp::vector)> rb_grad = [](pp::vector v) {
        pp::vector f_grad(2);
        double x = v[0]; double y = v[1];

        f_grad[0] = 2*(200*x*x*x - 200*x*y + x - 1);  //dfdx
        f_grad[1] = 200*(y-x*x); //dfdy

        return f_grad;
    };
    pp::vector rb_v(2); rb_v[0] = -2.0; rb_v[1] = 2.0;
    pp::vector root_rb_grad = pp::newton(rb_grad, rb_v);
    std::cout<<"Rosenbrock's valley function, extremum:\n\tx="<< root_rb_grad[0] <<"\n\ty="<<root_rb_grad[1]<<"\n";


    //Himmelblau's function
    // f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
    // df/dx = 2 * (2x * (x^2 + y - 11) + x + y^2 -7)
    // df/dy = 2 * (x^2 + 2y * (x + y^2 - 7) + y - 11)

    std::function<pp::vector(pp::vector)> hb_grad = [](pp::vector v) {
        pp::vector f_grad(2);
        double x = v[0]; double y = v[1];

        f_grad[0] = 2 * (2*x * (x*x + y - 11) + x + y*y -7);  //dfdx
        f_grad[1] = 2 * (x*x + 2*y * (x + y*y - 7) + y - 11); //dfdy

        return f_grad;
    };
    pp::vector hb_v(2); hb_v[0] = -2.0; hb_v[1] = 2.0;
    pp::vector root_hb_grad = pp::newton(hb_grad, hb_v);
    std::cout<<"Himmelblau's function, minimum:\n\tx="<< root_hb_grad[0] <<"\n\ty="<<root_hb_grad[1]<<\
               "\n\t(the function has four identical minima (of 0.0), I've just found one)\n";



    std::cout<<"\n\nEXERCISE B\n\n";
    //quick debugging to see if M works as intended
    auto [rlist1, ylist1] = M(-1, 1e-3, 8); double ME1 = ylist1.back()[0];
    auto [rlist2, ylist2] = M(-0.5, 1e-3, 8); double ME2 = ylist2.back()[0];
    auto [rlist3, ylist3] = M(0.0, 1e-3, 8); double ME3 = ylist3.back()[0];
    std::cout << "Debugging:\n";
    std::cout << "E = -1   M(E) = " << ME1 << "\n";
    std::cout << "E = -0.5 M(E) = " << ME2 << "\n";
    std::cout << "E =  0   M(E) = " << ME3 << "\n";


    pp::vector E0(1); E0[0] = -0.55;  // inital guess, should be -0.5
    double compE0 = compute_E(E0);
    std::cout<<"Finding the lowest root (E0): "<< compE0 <<"\n";

    // comparing wavefunction to exact value
    std::ofstream out("wavefunction.txt");
    auto [rlist, ylist] = M(compE0, 1e-3, 8);
    for(size_t i = 0; i < rlist.size(); i++){
        out << rlist[i] << " " << ylist[i][0] << " " << true_wavefunction(rlist[i]) << "\n";
    }

    //checking convergence
    std::cout << "\nChecking for convergence, E0=-0.5\n";
    std::cout << "Standard values: rmin=0.001, rmax=8, acc=0.01, eps=0.01\n";
    std::cout << "Altering rmin:\n";
    for(double rmin = 1; rmin > 2e-3; rmin/=5){
        double rminE = compute_E(E0, rmin, 8);
        std::cout << "\trmin=" << rmin << ", E=" << rminE << "\n";
    }
    std::cout << "Altering rmax:\n";
    for(int rmax = 2; rmax <= 8; rmax+=2){
        double rmaxE = compute_E(E0, 1e-3, rmax);
        std::cout << "\trmax=" << rmax << ", E=" << rmaxE << "\n";
    }
    std::cout << "Altering accuracy:\n";
    for(double acc = 1; acc > 2e-3; acc/=5){
        double accE = compute_E(E0, 1e-3, 8, acc);
        std::cout << "\tacc=" << acc << ", E=" << accE << "\n";
    }
    std::cout << "Altering epsilon:\n";
    for(double eps = 1; eps > 2e-3; eps/=5){
        double epsE = compute_E(E0, 1e-3, 8, 0.01, eps);
        std::cout << "\teps=" << eps << ", E=" << epsE << "\n";
    }




    return 0;
}
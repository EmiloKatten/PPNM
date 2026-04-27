#include<iostream>
#include<fstream>
#include"newton.h"


// Breit-Wigner function
double F(double E, double m, double Gamma, double A) {
    return A / ( (E - m)*(E - m) + (Gamma*Gamma)/4.0 );
};


int main(){

    std::cout<<"EXERCISE A\n\n";

    //Rosenbrock's valley function
    // f(x,y) = (1-x)^2 + 100(y-x^2)^2
    std::function<double(const pp::vector&)> f = [](const pp::vector& v) {
        double x = v[0];
        double y = v[1];
        return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
    };
    pp::vector f_x = {2.0, -2.0};
    auto [f_res, f_steps] = pp::newton(f, f_x);
    f_res.print("Minimum of Rosenbrock's vally (x,y):");
    std::cout<<"\t Correct minima: (1,1)\n";
    std::cout<<"\t # of steps: "<<f_steps<<"\n";

    //Himmelblau's function
    // f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
    std::function<double(const pp::vector&)> g = [](const pp::vector& v) {
        double x = v[0];
        double y = v[1];
        return (x*x + y - 11)*(x*x + y - 11) + (x + y*y - 7)*(x + y*y - 7);
    };
    pp::vector g_x = {5.0, 5.0};
    auto [g_res, g_steps] = pp::newton(g, g_x);
    g_res.print("A minimum of Himmelblau's vally (x,y):");
    std::cout<<"\t Correct minima: (3,2)\n";
    std::cout<<"\t # of steps: "<<g_steps<<"\n";



    std::cout<<"\n\nEXERCISE B\n\n";

    // Read and store the data
    std::vector<double> energy, signal, error;
    double x, y, z;
    while (std::cin >> x >> y >> z) {
        energy.push_back(x);
        signal.push_back(y);
        error.push_back(z);
    }

    // Deviation function
    std::function<double(const pp::vector&)> D = [&](const pp::vector& v) {
        double m = v[0];
        double gamma = v[1];
        double A = v[2];
        double sum = 0.0;

        for (size_t i = 0; i < energy.size(); i++){
            double diff = (F(energy[i], m, gamma, A) - signal[i]) / error[i];
            sum += diff * diff;
        };
        return sum;
    };

    // Initial guess
    pp::vector x0 = {125.0, 3.0, 10.0};

    auto [res, steps] = pp::newton(D, x0);
    double m = res[0]; double gamma = res[1]; double A = res[2];
    std::cout << "Optimal parameters:\n\tm = "<<m<<"\n\tΓ = "<<gamma<<"\n\tA = "<<A<<"\n";
    std::cout << "# of steps: "<<steps<<"\n";

    // Create fit
    std::ofstream out("fit.txt");
    for (double e=energy[0]; e<energy[energy.size()-1]; e+=0.2){
        out << e << " " << F(e, m, gamma, A) <<"\n";
    };


    return 0;
}
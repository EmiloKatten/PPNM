#include"hamiltonian.h"

pp::vector make_r(double rmax, double dr){
    int npoints = (int)(rmax/dr)-1;
    pp::vector r(npoints);
    for(int i=0;i<npoints;i++){r[i]=dr*(i+1);}
    return r;
}

pp::matrix make_H(double rmax, double dr){
    pp::vector r = make_r(rmax, dr);
    int npoints=r.size();
    pp::matrix H(npoints,npoints);
    for(int i=0;i<npoints-1;i++){
        H[i,i]  =-2*(-0.5/dr/dr);
        H[i,i+1]= 1*(-0.5/dr/dr);
        H[i+1,i]= 1*(-0.5/dr/dr);
    }
    H[npoints-1,npoints-1]=-2*(-0.5/dr/dr);
    for(int i=0;i<npoints;i++)H[i,i]+=-1/r[i];
    return H;
}

void hamiltonian(){
    double rmax = 30.0; 
    double dr = 0.1;
    pp::matrix H = make_H(rmax, dr);
    pp::vector r = make_r(rmax, dr);

    int n = H.size1();
    pp::vector f1(n+1);
    pp::vector f2(n+1);
    pp::vector f3(n+1);

    // The newton method was off by some factor, and this is due to different normalizations
    // With chatbot, I discovered the difference is sqrt(dr) (as also discussed in homework EVD)
    double factor = std::sqrt(dr);
    for(int i=0;i<n;i++){ // expected eigenfunctions and eigenvalues for n=1,2,3
        f1[i] = 2*r[i]*std::exp(-r[i])*factor;   
        f2[i] = 1/std::sqrt(2)*r[i]*(1-r[i]/2.0)
                * std::exp(-r[i]/2.0) * factor;
        f3[i] = 2.0/3.0/std::sqrt(3)*r[i]
                * (1-2*r[i]/3.0+2*r[i]*r[i]/27.0)
                * std::exp(-r[i]/3.0) * factor;
    }
    f1[n] = -1.0/2.0; // lambda n=1
    f2[n] = -1.0/8.0; // lambda n=2
    f3[n] = -1.0/18.0; // lambda n=3

    auto F = [&](pp::vector y){
        return eigen_f(H, y);
    };
    auto J = [&](pp::vector y){
        return eigen_J(H, y);
    };

    auto start_H = std::chrono::high_resolution_clock::now();
    auto [f1_newton, f1_steps] = newton(F, J, f1);
    auto [f2_newton, f2_steps] = newton(F, J, f2);
    auto [f3_newton, f3_steps] = newton(F, J, f3);
    auto stop_H = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_H = stop_H - start_H;

    std::ofstream out("wavefunctions.data");
    out << "0 0 0 0 0 0 0\n";
    for(int i=0;i<n;i++){
        out << r[i] << " " << f1_newton[i]*f1_newton[i] << " " << f2_newton[i]*f2_newton[i] << " " << f3_newton[i]*f3_newton[i] \
            << " " << f1[i]*f1[i] << " " << f2[i]*f2[i] << " " << f3[i]*f3[i] << "\n";
    }
    std::cout << "Using newton rootfinding, the eigenvalues of the hamiltonians (with rmax=30, dr=0.1):"\
              << "\n\tn=1 : " << f1_newton[n] << ",   error: "<< std::abs(f1_newton[n]-f1[n])\
              << "\n\tn=2 : " << f2_newton[n] << ",   error: "<< std::abs(f2_newton[n]-f2[n])\
              << "\n\tn=3 : " << f3_newton[n] << ",   error: "<< std::abs(f3_newton[n]-f3[n])<<"\n";
    std::cout << "Time for evaluation: " << elapsed_H << "s\n"; 
}

void instability(){
    double rmax = 30.0; 
    double dr = 0.1;
    pp::matrix H = make_H(rmax, dr);
    pp::vector r = make_r(rmax, dr);

    int n = H.size1();
    pp::vector f(n+1);
    pp::vector true_f(n+1);

    double factor = std::sqrt(dr);
    for (double eps=0.0; eps<=0.5; eps+=0.1){
        for(int i=0;i<n;i++){
            double f1 = 2*r[i]*std::exp(-r[i]) * factor;
            double f2 = 1/std::sqrt(2)*r[i]*(1-r[i]/2.0) * std::exp(-r[i]/2.0) * factor;
            double f3 = 2.0/3.0/std::sqrt(3)*r[i] * (1-2*r[i]/3.0+2*r[i]*r[i]/27.0) * std::exp(-r[i]/3.0) * factor;
            double noise = (2.0 * rand() / RAND_MAX - 1.0)/5; // [-0.1,0.1]

            f[i] = (1-eps)*f2 + eps * (f1 + f3)/2 + eps * noise;
            true_f[i] = f2;
        }

        double f1_lambda = -0.5;
        double f2_lambda = -0.125;
        double f3_lambda = -1.0/18.0;

        f[n] = (1-eps) * f2_lambda 
                + eps * (f1_lambda + f3_lambda)/2;

        auto F = [&](pp::vector y){
            return eigen_f(H, y);
        };
        auto J = [&](pp::vector y){
            return eigen_J(H, y);
        };
        auto [f_newton, f_steps] = newton(F, J, f);
        double error = f_newton[n]-f2_lambda;
        std::cout << "\teps="<<eps<<": E_newton="<< f_newton[n] << ", error: "<< std::abs(error) <<", E_guess: "<<f[n]<<". ncalls: " << f_steps << "\n";

        // plot last iteration for eps=0.5
        if (std::abs(eps - 0.5) < 1e-12){
            std::ofstream out("perturbation.data");
            out << "0 0 0\n";
            for(int i=0;i<n;i++){
            out << r[i] << " " << f_newton[i] << " " << f[i] <<" "<< true_f[i] << "\n";
            }
        }
    }
}
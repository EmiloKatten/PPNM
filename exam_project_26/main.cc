#include<random>
#include<iostream>
#include<fstream>
#include<chrono>
#include"root.h"
#include"run_timing.h"
#include"hamiltonian.h"

void run_timing();
void hamiltonian();

int main(){

    /*  
        First, quick debugging
        Himmelblau's function
        f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
        This becomes 0 exactly when both squares are zero, hence:
        F(x,y) = {x^2 + y - 11, x + y^2 - 7}
        and
        J = {{2x, 1}, {1, 2y}}
    */
    std::function<pp::vector(pp::vector)> hb_F = [](pp::vector v) {
        pp::vector F(2);

        double x = v[0]; 
        double y = v[1];

        F[0] = x*x + y - 11;
        F[1] = x + y*y - 7;

        return F;
    };
    std::function<pp::matrix(pp::vector)> hb_J = [](pp::vector v) {
        pp::matrix J(2,2);

        double x = v[0]; 
        double y = v[1];

        J(0,0) = 2*x;
        J(0,1) = 1;
        J(1,0) = 1;
        J(1,1) = 2*y;

        return J;
    };
    pp::vector hb_v0(2); 
    hb_v0[0] = -2.0; 
    hb_v0[1] = 2.0;

    auto [root_hb,steps_hb] = pp::newton(hb_F, hb_J, hb_v0);
    std::cout<<"===DEBUGGING===\n\n";
    std::cout<<"For debugging, let's find a root of the Himmelblau's function using analytical Jacobian:"\
               "\n\tx="<< root_hb[0] <<"\n\ty="<<root_hb[1]<<\
               "\n(the function has four identical minima (of 0.0), I've just found one)\n";
    



    /*
        Investigating the scaling of the updated rootfinding
    */
    std::cout<<"\n\n\n===TIMING===\n\n";
    std::cout << "Investigating the scaling of rootfinding using analytical Jacobian.\n\n";
    std::cout << "From the lecture notes, newton rootfinding scales as O(n^2), while QR-decomposition\n"\
                 <<"scales as O(n^3) for a symmetric matrix. Thus we should expect a scaling of O(n^3).\n\n";
    run_timing(false);
    std::cout << "\nRealistically, the size of the matrix should be no larger than N=300, since\n"\
              << "it becomes too expensive and time consuming.\n";

    

    /*
        Computing several lowest state eigenfunctions of the Hydrogen atom
    */
    std::cout<<"\n\n\n===HYDROGEN HAMILTONIAN===\n\n";
    std::cout << "Using the optimized newton rootfinding, let's try\n"\
              << "to compute several lower states of the Hydrogen atom.\n\n";
    std::cout << "For hydrogen, the energy in atomic units is: \n\tE_n = -1/(2n^2),\n"\
              << "where n is the shell number. Taking the three lowest states:\n"\
              << "\tE_1 = -1/2 = -0.5\n\tE_2 = -1/8 = -0.125\n\tE_3 = -1/18 = -0.0555\n\n";
    hamiltonian();
    

    
    /*
        Using an optimized alpha based on minimizing the norm in newton rootfinding wrt. alpha
    */
    std::cout<<"\n\n\n===OPTIMIZED LINE SEARCH===\n\n";
    std::cout << "Doing a taylor expansion of F(x+alpha*dx) around x yields F(x) + alpha*J(x)*dx\n"
              << "Minimizing the linearized prediction phi(alpha) = ||F(x) + alpha*J(x)*dx||^2\n"\
              << "provides an optimal alpha = - dot(F(x),J*dx) / dot(J*dx,J*dx) \n\n";
    std::cout << "Timing again using optimized alpha:\n";
    run_timing(true);

    return 0;
}
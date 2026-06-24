#include<iostream>
#include"root.h"
#include"run_timing.h"
#include"hamiltonian.h"

void run_timing();
void run_compare_timing();
void hamiltonian();
void instability();

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
    run_timing();
    std::cout << "\nRealistically, the size of the matrix should be no larger than N=300, since\n"\
              << "it becomes too expensive and time consuming.\n";
    std::cout << "See figure timing.svg for a plot of the time dependence.\n";

    
    /*
        Comparing newton and EVD
    */
    std::cout<<"\n\n\n===NEWTON VS EVD===\n\n";
    std::cout<< "Newton rootfinding only finds one eigenvector and eigenvalue, while\n"\
             << "requiring good initial guesses. On the contrary, EVD does not require\n"\
             << "any guesses and finds all eigenpairs very efficiently. Both should\n"
             << "scale as O(n^3). Comparing these directly:\n\n";
    run_compare_timing();
    std::cout<< "\nEVD is slower, but has the advantage of finding all eigenpairs without\n"\
             << "requiring guesses. If we normalized the EVD times by matrix size n,\n"\
             << "EVD is vastly superior to newton per eigenpair. However, for very large\n"
             << "matrices, if only one eigenpair is required and you have an inital guess,\n"\
             << "newton may be advantagous.\n";
    std::cout<< "comp_timing.svg plots the comparison.\n";



    /*
        Computing several lowest state eigenfunctions of the Hydrogen atom
    */
    std::cout<<"\n\n\n===HYDROGEN HAMILTONIAN===\n\n";
    std::cout << "Using the optimized newton rootfinding, let's try\n"\
              << "to compute several lower states of the Hydrogen atom.\n\n";
    std::cout << "For hydrogen, the energy in atomic units is: \n\tE_n = -1/(2n^2),\n"\
              << "where n is the shell number. Taking the three lowest states:\n"\
              << "\tE_1 = -1/2 = -0.5\n\tE_2 = -1/8 = -0.125\n\tE_3 = -1/18 = -0.0555\n"\
              << "These will be our initial eigenvalue guesses. The intital guesses for\n"\
              << "eigenvectors are the well-established wavefunctions for s1, s2 and\n"\
              << "s3 of hydrogen.\n\n";
    hamiltonian();
    std::cout << "See hamiltonian.svg for the probability density found with newton vs analytical solution.\n"; 


    
    /*
        Checking the instability of newton rootfinding.
    */
    std::cout<<"\n\n\n===NEWTON INSTABILITY===\n\n";
    std::cout << "Rootfinding with newton requires a good initial guess.\n"
              << "To test instability to perturbation, we'll try to perturbe the\n"\
              << "analytical guess for hydrogen s2, by introducing components of \n"\
              << "the eigenvectors for s1 and s3. Random noise will also be added.\n"\
              << "The initial guess of eigenvalue will also be perturbed by s1 and\n"\
              << "s3 eigenvalues.\n"\
              << "The resulting E_newton, error, E_guess, and number of calls:\n\n";
    instability();
    std::cout << "\nWe can clearly see, as epsilon increases, more steps are required\n"\
              << "for newton to converge. The eigenvalue guess increases, but newton\n"\
              << "still converges to the correct eigenvalue.\n"\
              << "Only for epsilon=0.5, the eigenpair guess makes the newton method fall\n"\
              << "into the wrong basin. The perturbation.svg figure shows this.\n";

    return 0;
}
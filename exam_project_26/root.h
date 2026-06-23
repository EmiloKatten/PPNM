#pragma once
#include"qr.h"
#include<functional>

namespace pp{

    /* New routine - with help from chatbot */
    // inline vector newton(std::function<vector(vector)> F, std::function<matrix(vector)> J_, vector x, double acc=1e-7){
    //     vector fx=F(x);
    //     while(fx.norm()>acc){
    //         matrix J=J_(x);
    //         qr QRJ(J);
    //         vector Dx=QRJ.solve(-fx);
    //         x += Dx;
    //         fx=F(x);
    //     }
    //     return x;
    // }
    /* New routine, using analytical Jacobian */
    inline std::pair<vector,int> newton(std::function<vector(vector)> F, std::function<matrix(vector)> J_,\
                         vector x, double acc=1e-7, double alpha_min=1e-3, double max_iter = 1e3){
        vector fx = F(x);
        int steps = 0;
        for (int i = 0; i<max_iter; i++){
            if (fx.norm() < acc){break;}
            matrix J = J_(x);
            qr QRJ(J);
            vector Dx = QRJ.solve(-fx);
            double alpha = 1.0;
            vector z, fz;
            while (true) {
                z = x+alpha*Dx;
                fz = F(z);
                if (fz.norm() < (1-alpha/2)*fx.norm()) {break;}
                if (alpha < alpha_min) {break;}
                alpha/=2;
            }
            x = z;
            fx = fz;
            steps += 1;
        }
        return {x,steps};
    };

    // The function F, using a symmetric matrix A and vector x = {v, lambda}
    inline vector eigen_f(const matrix& A, const vector& x){
        int n = A.size1();

        vector v(n);
        for (int i=0; i<n; i++){
            v[i] = x[i];
        }
        double lambda = x[n];

        vector F(n+1);
        vector Av = A*v;
        for (int i=0; i<n; i++){
            F[i] = Av[i] - lambda*v[i];
        }
        F[n] = dot(v,v) - 1.0;

        return F;
    };

    // The Jacobian of F, analytically determined
    inline matrix eigen_J(const matrix& A, const vector& x){
        int n = A.size1();

        vector v(n);
        for (int i=0; i<n; i++){
            v[i] = x[i];
        }
        double lambda = x[n];

        matrix J(n+1,n+1);
        // A - lambda*I
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                J(i,j) = A(i,j);
            }
        }
        for (int i=0; i<n; i++){
            J(i,i) -= lambda;
        }
        // -v
        for (int i=0; i<n; i++){
            J(i, n) -= v[i];
        }
        // 2v
        for (int j=0; j<n; j++){
            J(n, j) = 2*v[j];
        }
        // 0
        J(n,n) = 0;

        return J;
    };


    // Constructing a function to return the eigenvalues of a symmetric matrix A:
    inline std::tuple<double,vector,int> eigen_newton(const matrix& A, double lambda_start, vector& v_start){
        int n = A.size1();
        vector x(n+1);
        for (int i=0; i<n; i++){
            x[i] = v_start[i];
        }
        x[n] = lambda_start;

        // making lambda functions to compute F and J
        auto F = [&](vector y){
            return eigen_f(A, y);
        };
        auto J = [&](vector y){
            return eigen_J(A, y);
        };

        // use analytic newton rootfinding to update x:
        auto [x_res, steps] = newton(F, J, x);
        x = x_res;
        
        // return eigenvector and eigenvalue
        vector v(n);
        for (int i=0; i<n; i++){
            v[i] = x[i];
        }
        double lambda = x[n];

        return std::make_tuple(lambda, v, steps);
    }; 



    // Pick the best alpha that minimizes the linearized prediction of the residual norm - made with chatbot
    inline double opt_alpha(const vector& f, const vector& Jp, double alpha_min){
        double denom = pp::dot(Jp, Jp);
        if (denom == 0.0) return 1.0;

        double alpha = -(pp::dot(f, Jp)) / denom;

        // safeguard
        if (!std::isfinite(alpha)) alpha = 1.0;

        // clamp to reasonable range
        alpha = std::max(alpha, alpha_min);
        alpha = std::min(alpha, 1.0);

        return alpha;
    };
    inline std::pair<vector, int> opt_newton(std::function<vector(vector)> F, std::function<matrix(vector)> J_,\
                                             vector x, double acc=1e-7, double alpha_min=1e-3, double max_iter = 1e3){
        vector fx = F(x);
        int n_steps = 0;
        for (int i = 0; i<max_iter; i++){
            if (fx.norm() < acc){break;}
            matrix J = J_(x);
            qr QRJ(J);
            vector Dx = QRJ.solve(-fx);
            vector Jp = J * Dx;

            double alpha = opt_alpha(fx, Jp, alpha_min);

            vector z = x+alpha*Dx;
            vector fz = F(z);
            x = z;
            fx = fz;
            n_steps+=1;
        }
        return {x, n_steps};
    };
    inline std::tuple<double,vector,int> opt_eigen_newton(const matrix& A, double lambda_start, vector& v_start){
        int n = A.size1();
        vector x(n+1);
        for (int i=0; i<n; i++){
            x[i] = v_start[i];
        }
        x[n] = lambda_start;

        // making lambda functions to compute F and J
        auto F = [&](vector y){
            return eigen_f(A, y);
        };
        auto J = [&](vector y){
            return eigen_J(A, y);
        };

        // use analytic newton rootfinding to update x:
        auto [x_res,steps] = opt_newton(F, J, x);
        x = x_res;
        // return eigenvector and eigenvalue
        vector v(n);
        for (int i=0; i<n; i++){
            v[i] = x[i];
        }
        double lambda = x[n];

        return std::make_tuple(lambda, v, steps);
    }; 
}
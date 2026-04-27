#pragma once
#include<functional>
#include<cmath>
#include<tuple>
#include"qr.h"

namespace pp
{

    vector gradient(std::function<double(const vector&)> phi, vector x){
        double phi_x = (phi(x));
        vector g_phi(x.size());
        for (int i=0; i<x.size(); i++){
            double dx_i = (1+std::fabs(x[i]))*std::pow(2,-26);
            x[i]+=dx_i;
            g_phi[i] = (phi(x)-phi_x)/dx_i;
            x[i]-=dx_i;
        }
        return g_phi;
    };

    matrix hessian(std::function<double(const vector&)> phi, vector x){
        matrix H(x.size(),x.size());
        vector g_phi_x = gradient(phi, x);
        for (int j=0; j<x.size(); j++){
            double dxj = (1+std::fabs(x[j]))*std::pow(2,-13);
            x[j]+=dxj;
            vector dg_phi = gradient(phi, x) - g_phi_x;
            for (int i=0; i<x.size();i++){
                H[i,j]=dg_phi[i]/dxj;
            }
            x[j]-=dxj;
        };
        return H;
    };

    std::tuple<vector, int> newton(std::function<double(const vector&)> phi, vector x, double acc=1e-3){
        int n_steps = 0;
        while (true){
            n_steps += 1;
            vector g = gradient(phi, x);
            if (g.norm() < acc){break;};   //job done
            matrix H = hessian(phi, x);
            qr decomp(H);
            vector dx = decomp.solve(-g);
            double lambda = 1;
            while (lambda >= 1.0/1024.0) {
                if (phi(x + lambda * dx) < phi(x)){break;}
                lambda /= 2;
            }
            x += lambda * dx;
        };
        return {x, n_steps};
    };


} // namespace pp

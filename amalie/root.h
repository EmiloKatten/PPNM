#pragma once
#include"qr.h"
#include<cmath>
#include<functional>

namespace pp{

    matrix jacobian(std::function<vector(vector)> f, vector x, vector fx/* ,matrix& J */){
        int n = x.size();
        matrix J(n,n);
        for (int i = 0; i<n; i++){
            for (int j = 0; j<n; j++){
                double dxj=std::fabs(x[j])*std::pow(2,-26);
                x[j]+=dxj;
                vector df=f(x)-fx;
                J(i,j) = df[i]/dxj;
                x[j]-=dxj;
            }
        }
        return J;

    }

    vector newton(std::function<vector(vector)> f, vector x, double acc=1e-2, double alpha_min=1e-3, int max_iter=1e3){
        vector fx = f(x);
        for (int i = 0; i<max_iter; i++){
            if (fx.norm() < acc){break;}
            matrix J = jacobian(f, x, fx);
            qr QRJ(J);
            vector Dx = QRJ.solve(-fx);
            double alpha = 1.0;
            vector z, fz;
            while (true) {
                z = x+alpha*Dx;
                fz = f(z);
                if (fz.norm() < (1-alpha/2)*fx.norm()) {break;}
                if (alpha < alpha_min) {break;}
                alpha/=2;
            }
            x = z;
            fx = fz;
        }
        return x;
    };

}
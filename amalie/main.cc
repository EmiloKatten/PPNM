#include"root.h"


double secular(pp::vector d, pp::vector u, int p, double lambda){
    double f = d[p] + 2*u[p] - lambda;
    for (int k=0; k<n;k++){
        if (k==p){continue;}
        f -= u[k]*u[k]/(d[k]-lambda);
    }
    return f;
}

std::function<pp::vector(pp::vector)> make_f()

int main(){

    
    std::function<pp::vector(pp::vector)> debugf = [](pp::vector x) {
        double lambda = x[0];
        pp::vector y(1);
        y[0] = secular(d, u, p, lambda);
        return y;
    };

    
}
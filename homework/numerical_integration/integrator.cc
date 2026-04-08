#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include "integrator.h"

namespace pp{

double integrate(std::function<double(double)> f, double a, double b, double acc, double eps, double f2, double f3){
    double h = b-a;
    // accept infinite limits
    if (a == -inf && b == inf){
        auto g = [&](double t){
            double x = t/(1-t*t);
            return f(x)*(1+t*t)/(1-t*t)/(1-t*t);};
        return cc_integrate(g, -1, 1, acc, eps);
    }
    if (a == -inf && b != inf){
        auto g = [&](double t){
            double x = b + t/(1+t);
            return f(x) * 1/(1+t)/(1+t);
        };
        return cc_integrate(g, -1, 0, acc, eps);
    }
    if (b == inf && a != -inf){
        auto g = [&](double t){
            double x = a + t/(1-t);
            return f(x) * 1/(1-t)/(1-t);
        };
        return cc_integrate(g, 0, 1, acc, eps);
    }

    if (std::isnan(f2)){ f2 = f(a+2*h/6); f3 = f(a+4*h/6);}
    double f1 = f(a+h/6); double f4 = f(a+5*h/6);
    double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
    double q = (f1+f2+f3+f4)/4*(b-a);
    double err = std::abs(Q-q); double tol = acc + eps * std::abs(Q);
    if (err < tol){return Q;} 
    else return integrate(f, a, (a+b)/2, acc/std::sqrt(2), eps, f1, f2) + \
		        integrate(f, (a+b)/2, b, acc/std::sqrt(2), eps, f3, f4);    
}

double cc_integrate(std::function<double(double)> f, double a, double b, double acc, double eps){
    auto g = [&](double theta){
        double x = (a+b)/2 + (b-a)/2 * std::cos(theta);
        return f(x) * (b-a)/2 * std::sin(theta);
    };
    return integrate(g, 0, PI, acc, eps);
}

double erf(double z, double acc, double eps){
    if (z<0){return -erf(-z);}
    if (0<=z && z<=1){
        auto f = [](double x){return std::exp(-x*x);};
        return 2.0/std::sqrt(PI) * integrate(f,0,z, acc=acc, eps=eps);
    }
    else{ //(1<z)
        auto g = [&z](double t){return std::exp(-(z+(1-t)/t)*(z+(1-t)/t))/t/t;};
        return 1.0-2.0/std::sqrt(PI) * integrate(g,0,1, acc=acc, eps=eps);
    }   
}

bool approx(double a, double b, double acc, double eps){
	double diff = std::abs(a - b);
	if (diff <= acc) return true;  // absolute tolerance
	double max_ab = std::max(std::abs(a), std::abs(b));
	return diff <= eps * max_ab;   // relative tolerance
}

}

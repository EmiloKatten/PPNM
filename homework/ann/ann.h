#include<functional>
#include<cmath>
#include"gd.h"
#include"lcg.h"

namespace pp
{
    struct ann{
        int n; /* number of hidden neurons */
        std::function<double(double)> f; /* activation function */
        vector p; /* network parameters */

        // the ann uses gaussian wavelet as activation function
        ann(std::function<double(double)> f,int n) : f(f), n(n), p(3*n) {}
        ann(int n) : f([](double z) {return z*std::exp(-z*z);}), n(n), p(3*n) {
            lcg rng(4052026);
            for (int i=0; i<n; i++){
                p[i] = 2*rng.next() - 1;  //[-1, 1]
                p[n+i] = 0.1 + rng.next();
                p[2*n + i] = 2*rng.next() - 1;  //[-1, 1]
                p[0] = -1; p[1] = 1; // force edges
            }
        }

        double response(double x, const vector& p){
            double sum = 0;
            for (int i = 0; i<n; i++) {
                double a = p[i]; double b = p[n+i]; double w = p[2*n+i];
                double z = (x-a)/b;
                sum += f(z) * w;
            }
            return sum;
        }
        void train(const vector& x, const vector& y){
            int ncalls = 0;
            std::function<double(const vector&)> cost = [&](const vector& v) {
                ncalls += 1;
                double sum = 0.0;
                for (size_t i = 0; i < x.size(); i++){
                    double diff = response(x[i], v) - y[i];
                    sum += diff * diff;
                };
                return sum;
            };
            //p=amoeba::go(cost,p);
            p = gd::go(cost, p, 1e-4, 0.001, 10000);
            std::cerr<<ncalls;
        }

        double diff1(double x, const vector& p){
            double sum = 0;
            for (int i = 0; i<n; i++) {
                double a = p[i]; double b = p[n+i]; double w = p[2*n+i];
                double z = (x-a)/b;
                double dfdz = std::exp(-z*z) * (1 - 2*z*z);
                double dzdx = 1/b;
                double dfdx = dfdz * dzdx;
                sum += dfdx * w;
            }
            return sum;
        }
        double diff2(double x, const vector& p){
            double sum = 0;
            for (int i = 0; i<n; i++) {
                double a = p[i]; double b = p[n+i]; double w = p[2*n+i];
                double z = (x-a)/b;
                double dfdz2 = 2 * z * std::exp(-z*z) * (2*z*z - 3);
                double dzdx2 = 1/(b*b);
                double dfdx2 = dfdz2 * dzdx2;
                sum += dfdx2 * w;
            }
            return sum;
        }  
        double integral(double x, const vector& p){
            double sum = 0;
            for (int i = 0; i<n; i++) {
                double a = p[i]; double b = p[n+i]; double w = p[2*n+i];
                double z = (x-a)/b;
                sum += w * b * (-0.5 * std::exp(-z*z));     // b comes from dx = b * dz
            }
            return sum;
        }
        
    };

} // namespace pp
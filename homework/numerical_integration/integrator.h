#pragma once
#include <functional>
#include <cmath>
#include <limits>

namespace pp{

double const PI = std::numbers::pi;
double const inf = std::numeric_limits<double>::infinity();

double cc_integrate(std::function<double(double)> f, double a, double b, double acc=0.001, double eps=0.001); 
double integrate(std::function<double(double)> f, double a, double b, double acc=0.001, double eps=0.001, double f2=std::nan(""), double f3=std::nan(""));
double erf(double z, double acc=0.001, double eps=0.001);
bool approx(double a, double b, double acc=0.001, double eps=0.001);

}

#pragma once
#include<vector>

namespace pp {

int binsearch(const std::vector<double>& x, double z);
double linterp(const std::vector<double>& x, const std::vector<double>& y, double z);
double linterpInteg(const std::vector<double>& x, const std::vector<double>& y, double z);

struct qspline{
    const int n;
    std::vector<double> x,y,b,c;
    // constructor decleration
	qspline(const std::vector<double>& x, const std::vector<double>& y);
	double eval(double z);
	double deriv(double z);
	double integ(double z);
};

}
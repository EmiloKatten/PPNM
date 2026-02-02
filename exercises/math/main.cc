#include<iostream>
#include<cmath>
#include<numbers>
#include"sfuns.h"
int main(){
	
	double sqrt2 = std::sqrt(2.0);
	std::cout << "sqrt(2) = " << sqrt2 << "\n";
	
	double pow = std::pow(2, 1.0/5.0);
	std::cout << "2^(1/5) = " << pow << "\n";

	double pi = std::numbers::pi;
	double exp = std::exp(pi);
	std::cout << "exp(pi) = " << exp << "\n";

	double e = std::numbers::e;
	double pi_e = std::pow(pi, e);
	std::cout << "pi^e = " << pi_e << "\n";

	for(double x=1;x<=9;x+=1)
		std::cout << "fgamma("<<x<<")=" << sfuns::fgamma(x) << "\t" <<
		"lngamma("<<x<<")=" << sfuns::lngamma(x) << "\n";
	return 0;

}

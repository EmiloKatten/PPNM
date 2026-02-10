#include<iostream>
#include<complex>
#include<cmath>
#include<numbers>
#include"sfuns.h"

using complex=std::complex<double>;
//constexpr double  pi = 3.14159265358979324;
//constexpr double  e = 2.71828182845904523;
constexpr complex I = complex(0,1);

constexpr double pi = std::numbers::pi;
constexpr double e = std::numbers::e;

int main(){
	
	std::cout << "sqrt(2) = " << std::sqrt(2.0) << "\n";
	std::cout << "2^(1/5) = " << std::pow(2, 1.0/5.0) << "\n";
	std::cout << "exp(pi) = " << std::exp(pi) << "\n";
	std::cout << "pi^e = " << std::pow(pi,e) << "\n";

	std::cout << "log(i)=" << std::log(I)   <<"\n";
	std::cout << "   i^i=" << std::pow(I,I) <<"\n";
	std::cout << "   π^i=" << std::pow(pi,I) <<"\n";
	std::cout << "   e^i=" << std::pow(e,I) <<"\n";

	for(double x=1;x<=9;x+=1)
		std::cout << "fgamma("<<x<<")=" << sfuns::fgamma(x) << "\t" <<
		"lngamma("<<x<<")=" << sfuns::lngamma(x) << "\n";
	return 0;

}

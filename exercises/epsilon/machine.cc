#include<limits>
#include<cmath>
#include<iostream>
#include"machine.h"

void machine(void){
	float       f=1.0f; while((float)      (1.0f+f) != 1.0f){f/=2.0f;} f*=2.0f;
	double      d=1.0d; while((double)     (1.0d+d) != 1.0d){d/=2.0d;} d*=2.0d;
	long double l=1.0L; while((long double)(1.0L+l) != 1.0L){l/=2.0L;} l*=2.0L;
	std::printf("      float eps=%g\n",f);
	std::printf("     double eps=%g\n",d);
	std::printf("long double eps=%Lg\n",l);

	std::cout << std::numeric_limits<float>::epsilon() << "\n";
	std::cout << std::numeric_limits<double>::epsilon() << "\n";
	std::cout << std::numeric_limits<long double>::epsilon() << "\n";

	std::cout << "float check: " << std::pow(2, -23) << "\n";
	std::cout << "double check: " << std::pow(2, -52) << "\n";
}

#include<iostream>
#include<cmath>
#include<iomanip>
#include"compare.h"

bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9){
	double diff = std::abs(a - b);
        if (diff <= acc) return true;  // absolute tolerance
	double max_ab = std::max(std::abs(a), std::abs(b));
	return diff <= eps * max_ab;   // relative tolerance
}

void compare(){
double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
double d2 = 8*0.1;

std::cout << "d1==d2? " << (d1==d2 ? "true":"false") << "\n"; 

std::cout << std::fixed << std::setprecision(17);
std::cout << "d1=" << d1 << "\n";
std::cout << "d2=" << d2 << "\n";

std::cout << "Are they 'equal' (1: true, 0: false)? \t" << approx(d1, d2) << "\n";
}


#include"vec.h"
#include<iostream>

int main(){
	pp::vec<double> a {1, 2, 3};
	a.x = 6;
	std::cout << a.x << " " << a.y << " " << a.z << "\n";
	pp::vec<double> b {2, 3, 4};		
	std::cout << "Dot product of " << a << " and " << b << " is: " << pp::dot(a,b) << "\n";

	std::cout << "Cross product of "<< a << " and " << b << " is: " <<pp::cross(a,b) << "\n";

	pp::vec<float> c;
	std::cout << "Empty float type: " << c << "\n";
	
	pp::vec<double> alpha(a);
	std::cout << "Copy of first vector: " << alpha << "\n";

	pp::vec<double> beta(a+b);
	std::cout << "Copy of sum of vectors: " << beta << "\n";

	std::cout << a << " - " << b << " = " << a-b << "\n";
	
	std::cout << a << " * 2.31: a = " << (a * 2.31) << "\n";

	std::cout << a << " /= 2 = " << (a/=2) << "\n";
	
	std::cout << "Norm of new vector: ||"<< a <<"|| = " << pp::norm(a) << "\n";

	pp::vec<double> gamma {1, 2, 3};
	std::cout << "\n" << gamma << "\n";
	pp::vec<double> gamma1(gamma*(1+1e-10));
	std::cout << gamma << " + 1e-10 = " << gamma1 << "\n";
	std::cout << "Are they approximately the same? " << pp::approx(gamma, gamma1) << "\n";

    pp::vec<double> gamma2(gamma*(1+1e-5));
    std::cout << gamma << " + 1e-5 = " << gamma2 << "\n";
    std::cout << "Are they approximately the same? " << pp::approx(gamma, gamma2) << "\n";
    return 0;
}

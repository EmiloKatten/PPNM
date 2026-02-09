#include"vec.h"
#include<iostream>

int main(){
	pp::vec a {1, 2, 3};
	a.x = 6;
	std::cout << a.x << " " << a.y << " " << a.z << "\n";
	pp::vec b {2, 3, 4};		
	std::cout << "Dot product of " << a << " and " << b << " is: " << pp::dot(a,b) << "\n";
	std::cout << "Cross product of "<< a << " and " << b << " is: " <<pp::cross(a,b) << "\n";

	pp::vec c;
	std::cout << c << "\n";
	
	pp::vec alpha(a);
	std::cout << alpha << "\n";

	pp::vec beta(a+b);
	std::cout << beta << "\n";

	std::cout << a << "- " << b << "= " << a-b << "\n";
	
	std::cout << a << "* 2.31: a = " << (a * 2.31) << "\n";

	std::cout << a << "/= 2 = " << (a/=2) << "\n";
	
	std::cout << "Norm of new vector: ||"<< a <<"|| = " << pp::norm(a) << "\n";

	
	pp::vec gamma {1, 2, 3};
	std::cout << "\n" << gamma << "\n";
	pp::vec gamma1(gamma*(1+1e-10));
	std::cout << gamma << "+ 1e-10 = " << gamma1 << "\n";
	std::cout << "Are they approximately the same? " << pp::approx(gamma, gamma1) << "\n";

        pp::vec gamma2(gamma*(1+1e-5));
        std::cout << gamma << "+ 1e-5 = " << gamma2 << "\n";
        std::cout << "Are they approximately the same? " << pp::approx(gamma, gamma2) << "\n";	
	return 0;
}

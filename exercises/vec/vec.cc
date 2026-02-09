#include<iostream>
#include<cmath>
#include"vec.h"
#define SELF (*this)

namespace pp {

vec& vec::operator+=(const vec& other){
	x+=other.x;
	y+=other.y;
	z+=other.z;
	return (*this);
	}

vec& vec::operator-=(const vec& other){
	x-=other.x;
	y-=other.y;
	z-=other.z;
	return (*this);
	}

vec& vec::operator/=(double n){
	x/=n;
	y/=n;
	z/=n;
	return (*this);
	}

vec& vec::operator*=(double n){
	x*=n;
	y*=n;
	z*=n;
	return (*this);
	}

vec operator+(vec a, const vec& b){
	a += b;
	return a;
	}

vec operator-(vec a, const vec& b){
	a-= b;
	return a;
	}

vec operator*(const vec& a, double n){
	vec r=a;
	r*=n;
	return r;
	}

vec operator/(const vec& a, double n){
	vec r=a;
	r/=n;
	return r;
	}

double dot(const vec& a, const vec& b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
	}

vec cross(const vec& a, const vec& b){
	return vec(
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
		);
	}

double norm(const vec& a){
	return std::sqrt(dot(a, a));
	}

void vec::print(const std::string& s) const {
	std::cout << s << x << " " << y << " " << z << std::endl;
	}

bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9){
	double diff = std::abs(a - b);
	if (diff <= acc) return true;  // absolute tolerance
	double max_ab = std::max(std::abs(a), std::abs(b));
	return diff <= eps * max_ab;   // relative tolerance
}

bool approx(const vec& a, const vec& b){
	if(!approx(a.x,b.x))return false;
	if(!approx(a.y,b.y))return false;
	if(!approx(a.z,b.z))return false;
	return true;
	}

std::ostream& operator<<(std::ostream& os, const vec& v){
	os << "{ " << v.x << ", " << v.y << ", " << v.z << " } ";
        return os;
	}

}

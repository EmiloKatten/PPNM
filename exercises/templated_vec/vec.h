#pragma once
#include<string>
#include<iostream>
#include<cstdio>
#include<cmath>
namespace pp{
template<typename T>
struct vec{
	T x, y, z;

	// ctors
	vec(T a, T b, T c){ // parm. ctor
		printf("parametrized constructor called... \n");
		x=a; y=b; z=c;
		}
	vec() : vec(0,0,0) { // default ctor	
		printf("default constructor called... \n");
		}
	vec(const vec&)=default; // copy ctor: vec a=b;
	vec(vec&&)=default; // move ctor: vec a=b+c;
	
	// dtor
	~vec(){ printf("destructor called... \n");}

	// assignments
	vec& operator=(const vec&); // copy: a+b;
	vec& operator=(vec&&); // move: a=b+c

	// member operators
    vec& operator+=(const vec& other){
	    x+=other.x;
	    y+=other.y;
	    z+=other.z;
	    return (*this);
	}
	vec& operator-=(const vec& other){
        x-=other.x;
	    y-=other.y;
	    z-=other.z;
	    return (*this);
	}
	vec& operator*=(T n){
        x*=n;
	    y*=n;
	    z*=n;
	    return (*this);
    }
	vec& operator/=(T n){
        x/=n;
	    y/=n;
	    z/=n;
	    return (*this);
	}

	void print(const std::string& s="") const;

	// stream output
	//friend std::ostream& operator<<(std::ostream&, const vec&);
	};

// ostream operator<<
template <typename T>
std::ostream& operator<<(std::ostream& os, const vec<T>& v) {
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os;
	}

// non-member operators
template <typename T>
vec<T> operator+(vec<T> a, const vec<T>& b){
	a += b;
	return a;
	}
template <typename T>
vec<T> operator-(vec<T> a, const vec<T>& b){
	a-= b;
	return a;
	}
template <typename T>
vec<T> operator*(const vec<T>& a, T n){
	vec r=a;
	r*=n;
	return r;
	}
template <typename T>
vec<T> operator/(const vec<T>& a, T n){
	vec r=a;
	r/=n;
	return r;
	}

//dot product
template<typename T>
T dot(const vec<T>& a, const vec<T>& b){
        return a.x*b.x + a.y*b.y + a.z*b.z;
        }

//vector (cross) product
template <typename T>
vec<T> cross(const vec<T>& a, const vec<T>& b){
	return vec<T>(
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
		);
	}

//norm
template <typename T>
double norm(const vec<T>& a){
	return std::sqrt(dot(a, a));
	}

//approx
bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9);

template <typename T>
bool approx(const vec<T>& a, const vec<T>& b){
	if(!approx(a.x,b.x))return false;
	if(!approx(a.y,b.y))return false;
	if(!approx(a.z,b.z))return false;
	return true;
	}

//print for debugging
template <typename T>
void vec<T>::print(const std::string& s) const {
	std::cout << s << x << " " << y << " " << z << std::endl;
	}
}

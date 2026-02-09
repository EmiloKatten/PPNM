#pragma once
#include<string>
#include<cstdio>
namespace pp{
struct vec{
	double x,y,z;

	// ctors
	vec(double a, double b, double c){ // parm. ctor
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
	vec& operator+=(const vec&);
	vec& operator-=(const vec&);
	vec& operator*=(double);
	vec& operator/=(double);

	void print(const std::string& s="") const;

	// stream output
	friend std::ostream& operator<<(std::ostream&, const vec&);
	};
// non-member operators
vec operator-(vec, const vec&);
vec operator+(vec, const vec&);
vec operator*(const vec&, double);
vec operator/(const vec&, double);

//dot product
double dot(const vec& a, const vec& b);

//vector (cross) product
vec cross(const vec& a, const vec& b);

//norm
double norm(const vec& a);

//approx
bool approx(const vec& a, const vec& b);
}

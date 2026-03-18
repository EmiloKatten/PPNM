#include"qspline.h"
#include<vector>
#include<cassert>
#include<cmath>
#include<fstream>

using vec = std::vector<double>;

namespace pp
{
    
int binsearch(const vec& x, double z)
	{/* locates the interval for z by bisection */ 
	assert( z>=x[0] && z<=x[x.size()-1] );
	int i=0, j=x.size()-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
	}

double linterp(const vec& x, const vec& y, double z){
	int i=binsearch(x,z);
  	double dx=x[i+1]-x[i];
  	assert(dx>0);
  	double dy=y[i+1]-y[i];
  	return y[i]+dy/dx*(z-x[i]);
  	}

double linterpInteg(const vec& x, const vec& y, double z){
    int i = binsearch(x,z);
    double sum = 0.0;

    // full intervals
    for (int k=0; k<i; k++){
        sum += 0.5 * (x[k+1]-x[k])*(y[k]+y[k+1]);
    }

    // partial intervals
    double dx = x[i+1]-x[i];
    double dy = y[i+1]-y[i];
    double yz = y[i] + dy/dx*(z-x[i]);
    sum += 0.5*(z-x[i])*(y[i]+yz);
    return sum;
}

qspline::qspline(const vec& xx, const vec& yy) : x(xx), y(yy), n(xx.size()), b(n-1), c(n-1) {
    vec p(n-1), h(n-1);
		for(int i=0;i<n-1;i++){
			h[i]=x[i+1]-x[i];
			p[i]=(y[i+1]-y[i])/h[i];
		}
    //choose starting point and make recursion down:
	c[0]=0;
	for(int i=0;i<n-2;i++)
		c[i+1]=(p[i+1]-p[i]-c[i]*h[i])/h[i+1];
    //recursion down 
	c[n-2]/=2;
	for(int i=n-3;i>=0;i--)
		c[i]=(p[i+1]-p[i]-c[i+1]*h[i+1])/h[i];
    // determine b
	for(int i=0;i<n-1;i++)
		b[i]=p[i]-c[i]*h[i];
}
double qspline::eval(double z){
    assert(z >= x[0] && z <= x[x.size()-1]);
    int i = 0, j = n-1;
    while (j-i>1){
        int m = (i+j)*0.5;
        if (z > x[m]) i = m; else j = m;
    }
    double h = z - x[i];
    return y[i] + h*(b[i] + h*c[i]);
}
double qspline::deriv(double z){
    int i = binsearch(x, z);
    return b[i] + 2*c[i]*(z-x[i]);
}
double qspline::integ(double z){
    int i = binsearch(x,z);
    double sum = 0.0;

    // full intervals
    for (int k=0; k<i; k++){
        double h = x[k+1]-x[k];
        sum += y[k]*h + 0.5*b[k]*h*h + 1/3.0*c[k]*h*h*h;
    }

    // partial intervals
    double dx = z - x[i];
    sum += y[i]*dx + 0.5*b[i]*dx*dx + 1/3.0*c[i]*dx*dx*dx;
    return sum;
}

}
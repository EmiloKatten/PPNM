#include"matrix.h"

int main(){
    int n = 5;
    pp::vector v(n);
    pp::vector u(n);
    for(int i=0;i<n;i++)v[i]=i+1;
    for(int i=0;i<n;i++)u[i]=i+100;
    v.print("v=");
    u.print("u=");
    v+=u;
    v.print("new v=");
    u*=3;
    u.print("new u=");
    return 0;
}
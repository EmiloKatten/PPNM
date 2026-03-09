#pragma once
#include"matrix.h"
#include"qr.h"

namespace pp
{     
vector lsfit(const std::vector<std::function<double(double)>>& fs,
                const vector& x,
                const vector& y,
                const vector& dy){
    int n = x.size();
    int m = fs.size();

    matrix A(n,m);
    vector b(n);

    for(int i=0;i<n;i++){
        b[i] = y[i]/dy[i];

        for(int k=0;k<m;k++){
            A(i,k) = fs[k](x[i])/dy[i];
        }
    }

    qr qr(A);
    return qr.solve(b);
}
} // namespace pp

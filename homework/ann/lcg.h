#pragma once
#include <cmath>

namespace pp
{
    struct lcg{
        long seed; long a; long c; long m;
        double next(){
            seed = (a * seed + c) % m;
            return (double)seed/m;
        }
        lcg(long seed, long a=1664525, long c=1013904223, long m=std::pow(2,32)):
        seed(seed),a(a),c(c),m(m){};
    };
} // namespace pp

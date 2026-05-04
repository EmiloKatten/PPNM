#include"matrix.h"

/* gradient descent made with help from chat bot */

namespace pp {

struct gd {

static vector go(
    std::function<double(const vector&)> F,
    vector x,
    double acc = 1e-4,
    double step = 0.01,
    int maxsteps = 10000
){
    double eps = 1e-6;

    auto grad = [&](const vector& x){
        vector g(x.size());
        double fx = F(x);

        for (int i = 0; i < x.size(); i++){
            vector x1 = x;
            x1[i] += eps;
            g[i] = (F(x1) - fx) / eps;
        }
        return g;
    };

    for (int k = 0; k < maxsteps; k++){
        vector g = grad(x);

        if (g.norm() < acc) break;

        for (int i = 0; i < x.size(); i++){
            x[i] -= step * g[i];
        }
    }

    return x;
}

};

}
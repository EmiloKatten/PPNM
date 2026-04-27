#pragma once
#include"matrix.h"

namespace pp {

struct qr {
    matrix Q;
    matrix R;

    qr(const matrix& A)
        : Q(A), R(A.size2(), A.size2())
    {
        int m = A.size2();

        for (int i = 0; i < m; i++) {
            double norm = Q[i].norm();
            R(i, i) = norm;

            Q[i] /= norm;

            for (int j = i + 1; j < m; j++) {
                R(i, j) = pp::dot(Q[i], Q[j]);
                Q[j] -= Q[i] * R(i, j);
            }
        }
    }
	vector solve(vector b){
		int m = R.size1();  // R is m x m
		vector y = Q.transpose() * b;  // y = Q^T * b
		vector x(m);

		// Back substitution
		for (int i = m - 1; i >= 0; i--) {
			x[i] = y[i];
			for (int j = i + 1; j < m; j++) {
				x[i] -= R(i, j) * x[j];
			}
			x[i] /= R(i, i);
		}
		return x;
	}
	double det(){
		int m = R.size1();
		vector r(m);
		for (int i=0; i<m;i++){
			for (int j=0;j<m;j++){
				if (i == j){
					r[i] = R(i,i);
				}
			}
		}
		double prod = 1;
		for (int i=0; i<m;i++){
			prod *= r[i];
		}
		return prod;
	}
	matrix inverse() {
		int m = R.size1();
		matrix inv(m, m);

		for (int i = 0; i < m; i++) {
			vector e(m);
			e[i] = 1.0;

			vector col = solve(e);

			for (int j = 0; j < m; j++)
				inv(j, i) = col[j];  // store as column i
		}

		return inv;
	}
};

}
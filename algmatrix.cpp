#include "algmatrix.h"

double robust(double V, double rms)
{
	double k0 = 1.5;
	double k1 = 3.0;
	double scale = 1.0;

	double v = fabs(V / rms);
	if (v > k0 && v <= k1) {
		scale = (k0 / v) * pow((k1 - v) / (k1 - k0), 2);
	}
	else if (v > k1) {
		scale = 0;
	}

	return scale;
}

VecXd wlsq_LU(IN MatXd& A, IN VecXd& b, IN VecXd& W)
{
    // 计算加权矩阵和向量
    MatXd AWA = A.transpose() * W.asDiagonal() * A;
    VecXd AWb = A.transpose() * W.asDiagonal() * b;
    VecXd x = VecXd::Zero(A.cols());

    // 使用LU分解来求解线性方程
    if (AWA.fullPivLu().isInvertible()) {
        x = AWA.lu().solve(AWb);
    }

    return x;
}
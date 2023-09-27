#include "algmatrix.h"

VecXd wlsq_LU(IN MatXd& A, IN VecXd& b, IN VecXd& W, OUT MatXd* Q)
{
    // 计算加权矩阵和向量
    MatXd AWA = A.transpose() * W.asDiagonal() * A;
    VecXd AWb = A.transpose() * W.asDiagonal() * b;
    VecXd x = VecXd::Zero(A.cols());

    // 使用LU分解来求解线性方程
    if (AWA.fullPivLu().isInvertible()) {
        x = AWA.lu().solve(AWb);
    }

    // 计算X的协方差矩阵
    if (Q && Q->rows() == A.cols() && Q->cols() == A.cols()) {
        Q->resize(A.cols(), A.cols());
        (*Q) = AWA.inverse();
    }

    return x;
}

VecXd sliceVecByRate(IN VecXd& V)
{
    double vmean, vstd;
    int rows = V.rows();
    int nbad = 0, ibeg = 0, iend = 0, idx = 0;
    vector<double> vec;

    for (int i = 0; i < rows; i++) { vec.push_back(V(i)); }
    calcMeanStd(vec, vmean, vstd);

    for (int i = 0; i < rows; i++) {
        if (fabs(vec[i] - vmean) > 3 * vstd) {
            nbad++;
        }
    }
    if (nbad == 0)    { return V; }
    if (nbad == rows) { return VecXd(0); }

    stable_sort(vec.begin(), vec.end());
    for (int i = 0; i < nbad; ++i) {
        ibeg = i;
        iend = rows - (nbad - ibeg);
        if (fabs(vec[ibeg]) < fabs(vec[iend])) {
            break;
        }
        else {
            ibeg++;
        }
    }

    VecXd ret = VecXd::Zero(rows - nbad);
    for (int i = ibeg; i < iend; i++) {
        ret(idx++) = vec[i];
    }

    return ret;
}
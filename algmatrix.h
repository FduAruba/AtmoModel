#ifndef ALORITHM_MATRIX_H
#define ALORITHM_MATRIX_H

#include "comminterface.h"

double robust(double V, double rms);
VecXd wlsq_LU(IN MatXd& A, IN VecXd& b, IN VecXd& W);

#endif // !ALORITHM_MATRIX_H


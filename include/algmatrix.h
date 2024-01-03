#ifndef ALORITHM_MATRIX_H
#define ALORITHM_MATRIX_H

#include "comminterface.h"

VecXd LSQ(IN MatXd& A, IN VecXd& b, IN VecXd& W, OUT MatXd* Q);
VecXd sliceVecByRate(IN VecXd& V);

#endif // !ALORITHM_MATRIX_H


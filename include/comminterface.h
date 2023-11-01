#ifndef COMM_INTERFACE_H
#define COMM_INTERFACE_H

#include <Dense>
#include <SVD>
#include <LU>
#include "const.h"
#include "commtime.h"
#include "commfun.h"

typedef Eigen::MatrixXd MatXd;
typedef Eigen::VectorXd VecXd;
typedef map<string, map<int, FILE*>> FileFps;

#endif // !COMM_INTERFACE_H


#ifndef _DBVN_H_
#define _DBVN_H_

#include <RcppArmadillo.h>
#include <math.h>

double dbvnorm(double y1, double y2, double R);
vecmat Dbvn(double y1, double y2, double R);

#endif /* _DBVN_H_ */

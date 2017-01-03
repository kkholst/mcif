#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <RcppArmadillo.h>
#include "vmat.h"
#include "gmat.h"
#include "DataPairs.h"

double logdF1(unsigned row, const unsigned &cause, const unsigned &indiv, const DataPairs &data, const gmat &sigma, vec u);
double F1(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data, const gmat &sigma, vec u); // Marginal
double F1(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data); // Marginal, full follow-up
double F1(unsigned row, unsigned cause, unsigned indiv, unsigned cond_cause, const DataPairs &data, const gmat &sigma, vec u); // Conditional

double logdF2(unsigned row, const irowvec &causes, const DataPairs &data, const gmat &sigma, vec u);
double F2(unsigned row, irowvec causes, const DataPairs &data, const gmat &sigma, vec u);

#endif /* _FUNCTIONS_H_ */


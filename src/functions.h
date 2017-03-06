#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "vmat.h"
#include "gmat.h"
#include "DataPairs.h"
#include "pn.h"
#include "dbvn.h"

double logdF1(unsigned row, const unsigned &cause, const unsigned &indiv, const DataPairs &data, const gmat &sigma, vec u);
rowvec dlogdF1du(unsigned row, const unsigned &cause, const unsigned &indiv, const DataPairs &data, const gmat &sigma, vec u);

double F1(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data, const gmat &sigma, vec u); // Marginal
rowvec dF1du(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data, const gmat &sigma, vec u);

double F1(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data); // Marginal, full follow-up
rowvec dF1du(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data);

double F1(unsigned row, unsigned cause, unsigned indiv, unsigned cond_cause, const DataPairs &data, const gmat &sigma, vec u); // Conditional
rowvec dF1du(unsigned row, unsigned cause, unsigned indiv, unsigned cond_cause, const DataPairs &data, const gmat &sigma, vec u);

double logdF2(unsigned row, const irowvec &causes, const DataPairs &data, const gmat &sigma, vec u);
rowvec dlogdF2du(unsigned row, const irowvec &causes, const DataPairs &data, const gmat &sigma, vec u);

double F2(unsigned row, irowvec causes, const DataPairs &data, const gmat &sigma, vec u);
rowvec dF2du(unsigned row, irowvec causes, const DataPairs &data, const gmat &sigma, vec u);


#endif /* _FUNCTIONS_H_ */


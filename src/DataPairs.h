#ifndef DATAPAIRS_H
#define DATAPAIRS_H

#include "vmat.h"
#include "gmat.h"

class DataPairs { // Nyt navn?
public:
  unsigned ncauses;
  unsigned ncol;
  unsigned nrow;

  imat causes; // Causes 1,...,k, censoring=0, <0: t=tau (contribution to likelihood is just pi_j)

  mat alpha;
  mat dalpha;
  mat beta;
  mat gamma;
  mat pi;
  mat dpidu;
  mat dlogpidu;

  // Constructors
  // Constructor checker at der er lige mange rækker i alle matricer...
  DataPairs(unsigned ncauses, imat causes, mat alpha, mat dalpha, mat beta, mat gamma);
  DataPairs(unsigned ncauses);

  // Member functions
  irowvec causes_get(int i) const ;
  int causesMarg_get(int i, int indiv) const;

  rowvec alpha_get(int i, irowvec bothcauses) const;
  double alphaMarg_get(int i, int cause, int indiv) const;

  rowvec dalpha_get(int i, irowvec bothcauses) const;
  double dalphaMarg_get(int i, int cause, int indiv) const;

  rowvec beta_get(int i, irowvec bothcauses) const;
  double betaMarg_get(int i, int cause, int indiv) const;

  rowvec gamma_get(int i, irowvec bothcauses) const;
  double gammaMarg_get(int i, int cause, int indiv) const;

  rowvec pi_get(int i, irowvec bothcauses) const;
  double piMarg_get(int i, int cause, int indiv) const;

  double dpiduMarg_get(int i, int cause, int indiv, int dcause) const;
  rowvec dpiduMarg_get(int i, int cause, int indiv) const;

  double dlogpiduMarg_get(int i, int cause, int indiv, int dcause) const;
  rowvec dlogpiduMarg_get(int i, int cause, int indiv) const;

  void pi_gen(int i, vec u);
  void dpidu_gen(int i, vec u);
  void dlogpidu_gen(int i, vec u);
};

#endif /* DATAPAIRS_H */

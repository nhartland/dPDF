// $Id: parametrisation.h 3196 2015-08-27 14:16:27Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include <vector>
#include "NNPDF/common.h"
#include <gsl/gsl_vector.h>

using NNPDF::real;

class NostateMLP
{
public:
  NostateMLP(std::vector<int> const& arch);
  ~NostateMLP();

  void Compute(const gsl_vector* par, real const& x, real* out) const;
  int const& GetNParameters() const {return fNParameters;};
protected:
  const std::vector<int> fArch;
  const int fNParameters;
  const int fNOutputs;
  mutable real* fOutput;
};



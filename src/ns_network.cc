// $Id: parametrisation.cc 3195 2015-08-27 10:29:43Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <cstdlib>
#include <cmath>
#include <numeric>

#include "ns_network.h"

using namespace NNPDF;


int get_nparam(std::vector<int> const& arch)
{
  int nparam = 0;
  for (int i=1; i<arch.size(); i++)
    nparam += arch[i]*arch[i-1] + arch[i];
  return nparam;
}

NostateMLP::NostateMLP(std::vector<int> const& arch):
fArch(arch),
fNParameters(get_nparam(arch)),
fNOutputs(std::accumulate(arch.begin(), arch.end(), 0)),
fOutput(new real[fNOutputs])
{
  if (fArch[0] != 2) 
  {
    std::cerr << "Hard-coded for 2 values in the first layer" <<std::endl;
    exit(-1);
  }
}

NostateMLP::~NostateMLP()
{
  delete[] fOutput;
}

void NostateMLP::Compute(const gsl_vector* par, real const& x, real* out) const
{
  fOutput[0] = x;
  fOutput[1] = log(x);

  int ipar = 0, lout = 0;
  int iout = fArch[0];
  for (int i=1; i<fArch.size(); i++)
  {
    for (int j=0; j<fArch[i]; j++)
    {
      real h=0;
      for (int k=0; k<fArch[i-1]; k++)
        h-= gsl_vector_get(par,ipar++)*fOutput[lout+k];
      h+=gsl_vector_get(par,ipar++); // Bias
      fOutput[iout++] =  tanh(h);
    }
    lout += fArch[i];
  }

  for (int i=0; i<fArch[fArch.size()-1]; i++)
    out[i] = fOutput[fNOutputs - i - 1 ];
}
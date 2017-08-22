// $Id: minimizer.h 1286 2013-10-28 11:54:20Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iomanip>  

#include "NNPDF/common.h"
#include "deuteronset.h"
using std::vector;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include <NNPDF/experiments.h>
#include <NNPDF/lhapdfset.h>
using NNPDF::Experiment;
using NNPDF::LHAPDFSet;

// *************************************************************************************

class CMAESParam
{
public:
  CMAESParam(size_t const& _n, size_t const& _lambda);
  const size_t lambda;
  const size_t mu;
  const size_t n;
  size_t eigenInterval;
  double expN;
  double mu_eff;
  double csigma;
  double dsigma;
  double cc;
  double c1;
  double cmu;
  std::vector<double> wgts;
};

/**
 *  \class CMAESMinimizer
 *  \brief CMA-ES minimiser
 */
 
class CMAESMinimizer
{
public:
  CMAESMinimizer(int const& n, int const& lambda, double const& sigma);
  ~CMAESMinimizer();

  void Iterate(LHAPDFSet*, DeuteronSet*, gsl_vector*, vector<Experiment> const&);
  void NormVect(gsl_vector*) const; //!< Normally distributed random vector

private:
  std::vector<gsl_vector*> Mutation(DeuteronSet* pdf, const gsl_vector* m) const;
  gsl_vector* Recombination(gsl_vector* m, vector<size_t> const& rank, std::vector<gsl_vector*> const& yvals) const;

  void CSA(gsl_vector const* yavg);
  void CMA(int const& NIte, vector<size_t> const& rank, std::vector<gsl_vector*> const& yvals, gsl_vector const* yavg);

  void ComputeEigensystem();
  std::vector<real> ComputeErf(LHAPDFSet*, DeuteronSet*, vector<Experiment> const&);
  
protected:
  const CMAESParam fCMAES;
  size_t fIte;
  double fSigma;
  gsl_vector *fpsigma, *fpc;
  gsl_matrix *fC, *fBD, *finvC;  
  gsl_eigen_symmv_workspace *fwrkspc;
}; 

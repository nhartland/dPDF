// CMA-ES minimisation code
// Nathan Hartland 2017

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <numeric>
#include <algorithm>

#include "cmaes.h"

#include <NNPDF/randomgenerator.h>
using NNPDF::RandomGenerator;

// ************************* CMA-ES MINIMIZER *****************************

// Initialises parameters for CMA-ES minimiser
CMAESParam::CMAESParam(size_t const& _n, size_t const& _lambda):
  lambda(_lambda),
  mu(floor(lambda/2.0)),
  n(_n),
  eigenInterval(0.0),
  expN(0),
  mu_eff(0),
  csigma(0),
  dsigma(0),
  cc(0),
  c1(0),
  cmu(0),
  wgts(lambda,0)
{
  // Set expN
  expN = sqrt(n)*(1.0-1.0/(4.0*n) + 1.0/(21.0*n*n) );

  // Initialise w prime vector
  vector<double> wpr(lambda, 0);
  for (int i=0; i<lambda; i++)
    wpr[i] = log( (lambda + 1.0) / 2.0) - log(i+1);

  // Calculate weight sums
  double psumwgt = 0.0;   double nsumwgt = 0.0;
  double psumwgtsqr = 0.0;   double nsumwgtsqr = 0.0;
  for (int i=0; i<lambda; i++)
    if (i < mu) {psumwgt += wpr[i]; psumwgtsqr += wpr[i]*wpr[i]; }
    else        {nsumwgt += wpr[i]; nsumwgtsqr += wpr[i]*wpr[i]; }

  mu_eff = psumwgt*psumwgt/psumwgtsqr;
  const double mu_eff_minus = nsumwgt*nsumwgt/nsumwgtsqr;

  const double alpha_cov = 2.0;
  const double cmupr = alpha_cov*(mu_eff - 2.0 + 1.0/mu_eff)/(pow(n+2.0,2) + alpha_cov*mu_eff/2.0);

  // Set constants
  csigma = (mu_eff + 2.0) / (n + mu_eff + 5.0);
  dsigma = 1.0 + 2.0*fmax(0,(sqrt((mu_eff - 1.0)/(n + 1.0)))-1.0) + csigma;
  cc = (4.0 + mu_eff/n) / (n + 4.0 + 2.0*mu_eff/n );
  c1 = alpha_cov / ( pow(n + 1.3, 2.0) + mu_eff );
  cmu = std::min(1.0 - c1, cmupr);

  double sumwgtpos = 0.0;
  double sumwgtneg = 0.0;
  for (int i=0; i<lambda; i++)
    if (wpr[i] > 0) sumwgtpos += wpr[i];
    else sumwgtneg += fabs(wpr[i]);

  const double alpha_mu_minus = 1.0 + c1/cmu;
  const double alpha_mueff_minus = 1.0 + (2*mu_eff_minus)/(mu_eff + 2.0);
  const double alpha_posdef_minus = (1.0-c1-cmu)/(n*cmu);
  const double alpha_min = fmin(alpha_mu_minus, fmin(alpha_mueff_minus, alpha_posdef_minus));

  // Eigensystem solution interval
  eigenInterval = (lambda/(c1+cmu)/n)/10.0;
  
  // ********************************** Normalising weights  ****************************************

  for (int i=0; i<lambda; i++)
    wgts[i] = wpr[i]*( wpr[i] > 0 ? 1.0/sumwgtpos:alpha_min/sumwgtneg);


  // Test weight sum normalisation
  const double sumtestpos = std::accumulate(wgts.begin(), wgts.begin()+mu, 0.0);
  const double sumtestneg = std::accumulate(wgts.begin()+mu, wgts.end(), 0.0);

  std::cout << "CMA-ES Minimiser parameters initialised:" <<std::endl;
  std::cout << "n: "<< n <<" lambda: " <<lambda <<" mu: "<<mu<<" e_int: " << eigenInterval<<std::endl;
  std::cout << "csigma: " <<csigma <<" dsigma: " <<dsigma <<" E|N|: " << expN <<std::endl;
  std::cout << "cc: " << cc << " c1: " <<c1 << " cmu: " <<cmu <<std::endl;
  std::cout << "sumWpos: "<< sumtestpos <<" sumWneg == -alpha_min: "<<sumtestneg<<" == "<<-alpha_min<<std::endl;
  std::cout << "-c1/cmu == sumW: "<<-c1/cmu<<"  ==  "<<sumtestpos + sumtestneg <<std::endl;
  std::cout << std::endl;
}

/**
 * @brief CMAESMinimizer is the CMA-ES minimizer
 * @param settings the global NNPDFSettings
 */
CMAESMinimizer::CMAESMinimizer(int const& n, int const& lambda, double const& sigma):
fCMAES(n, lambda),
fIte(0),
fSigma(sigma),
fpsigma(gsl_vector_calloc( n )),
fpc(    gsl_vector_calloc( n )),
fC(     gsl_matrix_calloc( n, n )),
fBD(    gsl_matrix_calloc( n, n )),
finvC(  gsl_matrix_calloc( n, n )),
fwrkspc(gsl_eigen_symmv_alloc( n ))
{
  gsl_matrix_set_identity(fC);
}

CMAESMinimizer::~CMAESMinimizer()
{
  gsl_vector_free(fpsigma);
  gsl_vector_free(fpc);
  gsl_matrix_free(fC);
  gsl_matrix_free(fBD);
  gsl_matrix_free(finvC);
  gsl_eigen_symmv_free(fwrkspc);
}

void CMAESMinimizer::ComputeEigensystem()
{
    // Initialise matrices
    gsl_matrix *B = gsl_matrix_calloc( fCMAES.n, fCMAES.n );
    gsl_matrix *D = gsl_matrix_calloc( fCMAES.n, fCMAES.n );
    gsl_matrix *invD = gsl_matrix_calloc( fCMAES.n, fCMAES.n );

    gsl_matrix_set_zero (fBD);
    gsl_matrix_set_zero (finvC);

    // Calculate the eigensystem
    gsl_matrix* C = gsl_matrix_calloc( fCMAES.n, fCMAES.n );
    gsl_vector* E = gsl_vector_calloc( fCMAES.n );
    gsl_matrix_memcpy (C, fC);
    gsl_eigen_symmv (C, E, B, fwrkspc);

    // Compute condition number
    double min, max;
    gsl_vector_minmax (E, &min, &max);
    const double K = max/min;
    std::cout << "CMA - ConditionNumber: " << K <<std::endl;
    // Initialise D, invD
    for (size_t i=0; i<fCMAES.n; i++)
    {
      gsl_matrix_set(D,i,i, sqrt(gsl_vector_get(E,i)));
      gsl_matrix_set(invD,i,i, 1.0/sqrt(gsl_vector_get(E,i)));
    }

    // Compute BD, Cinv, use C as a temporary
    gsl_matrix_set_zero (C);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, B, D, 0.0, fBD); // BD
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, invD, B, 0.0, C ); // D^-1 * B^T
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, B, C, 0.0, finvC ); // B * D^-1 * B^T

    gsl_matrix_free(C);
    gsl_matrix_free(B);
    gsl_matrix_free(D);
    gsl_matrix_free(invD);
    gsl_vector_free(E);
}

void CMAESMinimizer::Iterate(gsl_vector* m, const CostComputer* cost)
{
  // First setup the required matrices
  if (fIte++ % fCMAES.eigenInterval == 0 )
    ComputeEigensystem();

  // Setup and mutate PDF members
  vector<gsl_vector*> yvals;
  vector<gsl_vector*> xvals;
  Mutation(xvals, yvals, m);

  // Compute ERF and rank members
  const vector<NNPDF::real> erf = (*cost)(xvals);
  vector<NNPDF::real> erf_srt = erf; std::sort(erf_srt.begin(), erf_srt.end());
  vector<size_t> irank_map(fCMAES.lambda,0); // Weight-ordered map to members (index is i)
  for (int i=0; i<fCMAES.lambda; i++)
    irank_map[std::distance(erf_srt.begin(), std::find(erf_srt.begin(), erf_srt.end(), erf[i]))] = i;

  // Compute weighted shift and set new mean
  gsl_vector* yavg = Recombination(m, irank_map, yvals);

  // ********************************** Adaptation  ****************************************
  CSA(yavg); CMA(fIte + 1, irank_map, yvals, yavg );

  for (auto i : yvals ) gsl_vector_free(i);
  gsl_vector_free(yavg);

};

void CMAESMinimizer::Mutation(std::vector<gsl_vector*>& xv, std::vector<gsl_vector*>& yv, const gsl_vector* m) const
{
  gsl_vector* z = gsl_vector_calloc(fCMAES.n); 
  for (size_t i=0; i<fCMAES.lambda; i++)
  {
    gsl_vector* x = gsl_vector_calloc(fCMAES.n);
    gsl_vector* y = gsl_vector_calloc(fCMAES.n);
    gsl_vector_set_zero (z); NormVect(z);
    gsl_vector_set_zero (y); gsl_blas_dgemv (CblasNoTrans, 1.0, fBD, z, 1.0, y);
    gsl_vector_set_zero (x); gsl_vector_memcpy (x, m); gsl_blas_daxpy (fSigma, y, x);
    xv.push_back(x);
    yv.push_back(y);
  }
  gsl_vector_free(z);
}

gsl_vector* CMAESMinimizer::Recombination(gsl_vector* m, vector<size_t> const& irank_map, std::vector<gsl_vector*> const& yvals) const
{
  // Compute average step
  gsl_vector* yavg = gsl_vector_calloc(fCMAES.n);
  for (int i=0; i<fCMAES.mu; i++)
    gsl_blas_daxpy (fCMAES.wgts[i], yvals[irank_map[i]], yavg);

  gsl_vector *newm  = gsl_vector_calloc( fCMAES.n );
  gsl_vector_memcpy(newm, m); 
  gsl_blas_daxpy (fSigma, yavg, newm);

  // Set new mean - this is a bit silly, should just daxpy straight into it
  gsl_vector_memcpy(m, newm);
  gsl_vector_free(newm);
  return yavg;
}

// Cumulative step-size adaptation
void CMAESMinimizer::CSA( gsl_vector const* yavg )
{
  const double alpha = sqrt(fCMAES.csigma*(2.0 - fCMAES.csigma)*fCMAES.mu_eff ); // Coeff of matrix multiply
  const double beta = (1.0-fCMAES.csigma); // Coeff of sum
  gsl_blas_dgemv (CblasNoTrans, alpha, finvC, yavg, beta, fpsigma);
  double pnorm = 0; gsl_blas_ddot (fpsigma, fpsigma, &pnorm);

  const double sigrat = fCMAES.csigma/fCMAES.dsigma;
  fSigma = fSigma*exp(sigrat*(sqrt(pnorm)/fCMAES.expN - 1.0));

  std::cout << "CSA - StepSize: "<<fSigma <<" expFac: "<<sqrt(pnorm)/fCMAES.expN << std::endl;
}

// Covariance matrix adaptation
void CMAESMinimizer::CMA( int const& NIte, vector<size_t> const& irank_map, std::vector<gsl_vector*> const& yvals, gsl_vector const* yavg )
{
  // Compute norm of p-sigma
  const double pnorm = gsl_blas_dnrm2 (fpsigma);
  const double hl = pnorm / (sqrt(1.0 - pow(1.0 - fCMAES.csigma,2*(NIte+1))));
  const double hr = (1.4 + 2.0/(fCMAES.n + 1))*fCMAES.expN;
  const double hsig = (hl < hr) ? 1:0;
  const double dhsig = (1 - hsig)*fCMAES.cc*(2-fCMAES.cc);
  
  const double alpha = hsig*sqrt(fCMAES.cc*(2.0-fCMAES.cc)*fCMAES.mu_eff);
  gsl_vector_scale( fpc, (1.0-fCMAES.cc));
  gsl_blas_daxpy (alpha, yavg, fpc);

  const double weightsum = std::accumulate(fCMAES.wgts.begin(),fCMAES.wgts.end(), 0.0 );
  const double Cscale = (1.0 + fCMAES.c1*dhsig - fCMAES.c1 - fCMAES.cmu*weightsum );

  if ( Cscale != 1.0 ) gsl_matrix_scale(fC, Cscale);
  gsl_blas_dger (fCMAES.c1, fpc, fpc, fC); // Rank-1 update

  // Rank-mu update
  for (int i=0; i<fCMAES.lambda; i++)
  {
    const gsl_vector* yval = yvals[irank_map[i]];
    double wo = fCMAES.wgts[i];
    if (fCMAES.wgts[i] < 0)
    {
      gsl_vector *cy  = gsl_vector_calloc( fCMAES.n );
      gsl_blas_dgemv (CblasNoTrans, 1.0, finvC, yval, 1.0, cy);
      const double norm = gsl_blas_dnrm2 (cy);
      wo *= fCMAES.n / (norm*norm);
      gsl_vector_free(cy);
    }
    gsl_blas_dger (fCMAES.cmu*wo, yval, yval, fC); 
  }
}

void CMAESMinimizer::NormVect(gsl_vector* vec) const
{
    for (size_t i=0; i<vec->size; i++)
      gsl_vector_set(vec, i, RandomGenerator::GetRNG()->GetRandomGausDev(1));
}

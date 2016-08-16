#include <config.h>
#include "DSWTN.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <JRmath.h>
#include <cmath>

using std::vector;

// #define LAMBDA(par) (*par[0])
#define ALPHA(par) (*par[0])
#define THETA(par) (*par[1])
#define D(par) (*par[2])
#define V(par) (*par[3])


namespace jags {
namespace wald {

DSWTN::DSWTN() : ScalarDist("dswtn", 4, DIST_PROPORTION)
{}

bool DSWTN::checkParameterValue (vector<double const *> const &parameters) const
{
    return  (true);
}

double DSWTN::dswtn(double t, vector<double const *> const &parameters) const
{
  double lambda = 1; // LAMBDA(parameters); 
  double alpha = ALPHA(parameters); 
  double v = V(parameters); 
  double d = D(parameters);
  double theta = THETA(parameters);
    
  double w;
  
  /* w = alpha * sqrt( lambda / (2 * M_PI * pow(t, 3) * (lambda * t * v + 1)) ) *
      1 / pnorm(d / sqrt(v), 0, 1, 1, 0) *
      exp( - (lambda * pow(d * t - alpha, 2)) / 
           (2 * t * (lambda * t * v + 1)) ) *
      pnorm( (lambda * alpha * v + d) / 
             (sqrt(lambda * t * pow(v, 2) + v) ), 0, 1, 1, 0);
  */
  
  w = log(alpha) + 
      0.5 * ( log(lambda) - 
              log(2) - log(M_PI) - 3 * log(t-theta) - log(lambda * (t-theta) * v + 1)) -
      log(pnorm(d / sqrt(v), 0, 1, 1, 0) ) -
      (lambda * pow(d * (t-theta) - alpha,2)) / (2 * (t-theta) * (lambda * (t-theta) * v + 1)) +
      log(pnorm((lambda * alpha * v + d) / (sqrt(lambda * (t-theta) * pow(v,2) + v)), 0, 1, 1, 0 ));
  
  return w;
}

double DSWTN::logDensity(double x, PDFType type,
       vector<double const *> const &parameters,
       double const *lbound, double const *ubound) const 
{
    double d = 0;

    d = dswtn(x , parameters);
    
    return d;
}

double DSWTN::randomSample(vector<double const *> const &parameters, 
         double const *lbound, double const *ubound,
         RNG *rng) const
{
    return fabs( rng->uniform() );
}

double DSWTN::typicalValue(vector<double const *> const &parameters,
         double const *lbound, double const *ubound) const
{
    return 0.5;
}

bool DSWTN::isDiscreteValued(vector<bool> const &mask) const
{
    return true;
}

} // namespace wald
} // namespace jags

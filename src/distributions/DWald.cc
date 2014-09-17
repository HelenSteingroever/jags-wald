#include <config.h>
#include "DWald.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <JRmath.h>
#include <cmath>

using std::vector;

#define LAMBDA(par) (*par[0])
#define ALPHA(par) (*par[1])
#define V(par) (*par[2])
#define D(par) (*par[3])

namespace wald {

DWald::DWald() : ScalarDist("dwald_trunc", 4, DIST_PROPORTION)
{}

bool DWald::checkParameterValue (vector<double const *> const &parameters) const
{
    return  (true);
}

double DWald::dwald_trunc(double t, vector<double const *> const &parameters) const
{
  double lambda = LAMBDA(parameters); 
  double alpha = ALPHA(parameters); 
  double v = V(parameters); 
  double d = D(parameters);

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
              log(2) - log(M_PI) - 3 * log(t) - log(lambda * t * v + 1)) -
      log(pnorm(d / sqrt(v), 0, 1, 1, 0) ) -
      (lambda * pow(d * t - alpha,2)) / (2 * t * (lambda * t * v + 1)) +
      log(pnorm((lambda * alpha * v + d) / (sqrt(lambda * t * pow(v,2) + v)), 0, 1, 1, 0 ));
  
  return w;
}

double DWald::logDensity(double x, PDFType type,
       vector<double const *> const &parameters,
       double const *lbound, double const *ubound) const 
{
    double d = 0;

    d = dwald_trunc(x , parameters);
    
    return d;
}

double DWald::randomSample(vector<double const *> const &parameters, 
         double const *lbound, double const *ubound,
         RNG *rng) const
{
    return fabs( rng->uniform() );
}

double DWald::typicalValue(vector<double const *> const &parameters,
         double const *lbound, double const *ubound) const
{
    return 0.5;
}

bool DWald::isDiscreteValued(vector<bool> const &mask) const
{
    return true;
}

}

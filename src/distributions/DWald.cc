#include <config.h>
#include "DWald.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <JRmath.h>
#include <cmath>

#include <gsl/gsl_specfunc.h>

using std::vector;

#define ALPHA(par) (*par[0])
#define LAMBDA(par) (*par[1])
#define NU(par) (*par[2])

namespace wald {

DWald::DWald() : ScalarDist("dwald", 4, DIST_PROPORTION)
{}

bool DWald::checkParameterValue (vector<double const *> const &parameters) const
{
    return(true);
}

double DWald::dwald(double t, vector<double const *> const &parameters) const
{
  double alpha = ALPHA(parameters);
  double lambda = LAMBDA(parameters);
  double nu = NU(parameters);

  double logd;

  logd = log(alpha) + .5 * (log(lambda) - log(2)-log(M_PI) - 3.*log(t)) + ( -lambda * pow((nu*t-alpha),2) / (2.*t));

  return(logd);

}

double DWald::logDensity(double x, PDFType type,
       vector<double const *> const &parameters,
       double const *lbound, double const *ubound) const 
{
  double logd = 0;

  logd = dwald(x , parameters);

  return logd;
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

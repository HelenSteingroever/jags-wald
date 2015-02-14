#include <config.h>
#include "DSWald.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <JRmath.h>
#include <cmath>

using std::vector;

#define ALPHA(par) (*par[0])
#define NU(par) (*par[1])
#define THETA(par) (*par[2])

namespace wald {

DSWald::DSWald() : ScalarDist("dswald", 3, DIST_PROPORTION)
{}

bool DSWald::checkParameterValue (vector<double const *> const &parameters) const
{
    return(true);
}

double DSWald::dswald(double t, vector<double const *> const &parameters) const
{
  double alpha = ALPHA(parameters);
  double nu = NU(parameters);
  double theta = THETA(parameters);

  double logd;

  logd = log(alpha) +  (-.5) * (log(2) + log(M_PI) + 3*log(t-theta)) 
         + ( -pow((alpha-nu*(t-theta)),2) / (2*(t-theta)) );

  return(logd);
}

double DSWald::logDensity(double x, PDFType type,
       vector<double const *> const &parameters,
       double const *lbound, double const *ubound) const 
{
    double logd = 0;

    logd = dswald(x , parameters);

    return logd;
}

double DSWald::randomSample(vector<double const *> const &parameters, 
         double const *lbound, double const *ubound,
         RNG *rng) const
{
    return fabs( rng->uniform() );
}

double DSWald::typicalValue(vector<double const *> const &parameters,
         double const *lbound, double const *ubound) const
{
    return 0.5;
}

bool DSWald::isDiscreteValued(vector<bool> const &mask) const
{
    return true;
}

}

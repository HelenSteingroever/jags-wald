#include <config.h>
#include "PSWald.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <JRmath.h>
#include <cmath>

using std::vector;

#define ALPHA(par) (*par[0])
#define NU(par) (*par[1])
#define THETA(par) (*par[2])

namespace jags {
namespace wald {

PSWald::PSWald() : ScalarDist("pswald_upper", 3, DIST_PROPORTION)
{}

bool PSWald::checkParameterValue (vector<double const *> const &parameters) const
{
    return(true);
}

double PSWald::logDensity(double x, PDFType type,
       vector<double const *> const &parameters,
       double const *lbound, double const *ubound) const 
{
  double p = 0;

  if(x<=THETA(parameters)) 
    p = 0;
  else 
    p = pnorm((NU(parameters)*(x-THETA(parameters))-ALPHA(parameters)) /
              sqrt((x-THETA(parameters))), 0,1,1,0) + 
        exp(2*ALPHA(parameters)*NU(parameters)) *
            pnorm(-(NU(parameters)*(x-THETA(parameters))+ALPHA(parameters)) /
            sqrt((x-THETA(parameters))), 0,1,1,0);

  return log(1-p);
}

double PSWald::randomSample(vector<double const *> const &parameters, 
         double const *lbound, double const *ubound,
         RNG *rng) const
{
    return fabs( rng->uniform() );
}

double PSWald::typicalValue(vector<double const *> const &parameters,
         double const *lbound, double const *ubound) const
{
    return 0.5;
}

bool PSWald::isDiscreteValued(vector<bool> const &mask) const
{
    return true;
}

} // namespace wald
} // namespace jags

#include <config.h>
#include "DSWald.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <JRmath.h>
#include <cmath>

using std::vector;

#define ALPHA(par) (*par[0])
#define THETA(par) (*par[1])
#define NU(par) (*par[2]) // xi[t]

namespace jags {
namespace wald {

DSWald::DSWald() : ScalarDist("dswald", 3, DIST_PROPORTION)
{}

bool DSWald::checkParameterValue (vector<double const *> const &parameters) const
{
    return(true);
}

double DSWald::logDensity(double x, PDFType type,
       vector<double const *> const &parameters,
       double const *lbound, double const *ubound) const 
{
    double logd = 0;

    if(x<=THETA(parameters)) 
      logd = JAGS_NEGINF;
    else 
      logd = log(ALPHA(parameters)) +  (-.5) * (log(2) + log(M_PI) + 3*log(x-THETA(parameters))) +
             ( -pow((ALPHA(parameters)-NU(parameters)*(x-THETA(parameters))),2) / (2*(x-THETA(parameters))) );

    return logd;

    // code template for R provided by Robert Miller and adapted thereafter
    /*
    double d = 0;

    if(x<=THETA(parameters)) 
      d = 0;
    else 
      d = ALPHA(parameters) * sqrt(2*M_PI*pow((x-THETA(parameters)),3)) * 
          exp(-pow((ALPHA(parameters)-NU(parameters)*(x-THETA(parameters))),2) / (2*(x-THETA(parameters))));

    return log(d);
    */

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

} // namespace wald
} // namespace jags

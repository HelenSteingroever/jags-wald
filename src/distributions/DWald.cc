#include <config.h>
#include "DWald.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <JRmath.h>
#include <cmath>

using std::vector;

#define ALPHA(par) (*par[0])
#define LAMBDA(par) (*par[1])
#define NU(par) (*par[2])

namespace wald {

DWald::DWald() : ScalarDist("dwald", 3, DIST_PROPORTION)
{}

bool DWald::checkParameterValue (vector<double const *> const &parameters) const
{
    return(true);
}

double DWald::logDensity(double x, PDFType type,
       vector<double const *> const &parameters,
       double const *lbound, double const *ubound) const 
{
    double logd = 0;

    logd = log(ALPHA(parameters)) + .5 * (log(LAMBDA(parameters)) - 
           log(2)-log(M_PI) - 3.*log(x)) + ( -LAMBDA(parameters) * 
           pow((NU(parameters)*x-ALPHA(parameters)),2) / (2.*x));

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

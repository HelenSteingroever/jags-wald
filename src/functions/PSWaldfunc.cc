#include <config.h>
#include "PSWaldfunc.h"

#include <JRmath.h>
#include <cmath>

using std::vector;
using std::string;

#define T(par) (*par[0])
#define ALPHA(par) (*par[1])
#define NU(par) (*par[2])
#define THETA(par) (*par[3])

namespace jags {
namespace wald {

PSWaldfunc::PSWaldfunc()
  :ScalarFunction("pswald", 4)
{}

bool PSWaldfunc::checkParameterValue(vector<double const *> const &args) const
{
  return (true);
}

double PSWaldfunc::evaluate(vector<double const *> const &args) const
{
  double p;

  p = T(args) - THETA(args);
  p = pnorm((NU(args)*p-ALPHA(args))/sqrt(p), 0,1,1,0) + 
      exp(2*ALPHA(args)*NU(args))*pnorm(-(NU(args)*p+ALPHA(args))/sqrt(p), 0,1,1,0);

  return(p);
}


} // namespace wald
} // namespace jags

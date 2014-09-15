#include <config.h>
#include "LogBernFun.h"
#include <util/nainf.h>

#include <cmath>

using std::vector;
using std::string;

#define x(par) (*args[0])
#define prob(par) (*args[1])

namespace bernoulli {

LogBernFun::LogBernFun() :ScalarFunction("logbern", 2)
{}

bool LogBernFun::checkParameterValue(vector<double const *> const &args) const
{
  return (x(par) == 0.0 || x(par) == 1.0 
          && (prob(par) <= 1.0 && prob(par) >= 0.0));
}

double LogBernFun::evaluate(vector<double const *> const &args) const
{
  double d = 0;
  if (x(par) == 1)
    d = prob(par);
  else if (x(par) == 0)
    d = 1 - prob(par);
  
  return d == 0 ? JAGS_NEGINF : log(d);
}

}

#include <config.h>
#include "DSWaldfunc.h"

#include <JRmath.h>
#include <cmath>

using std::vector;
using std::string;

#define T(par) (*par[0])
#define ALPHA(par) (*par[1])
#define NU(par) (*par[2])
#define THETA(par) (*par[3])

namespace wald {

DSWaldfunc::DSWaldfunc()
  :ScalarFunction("dswald", 4)
{}

bool DSWaldfunc::checkParameterValue(vector<double const *> const &args) const
{
  return (true);
}

double DSWaldfunc::evaluate(vector<double const *> const &args) const
{
  double d;

  if(T(args)<=THETA(args)) 
    d = 0;
  else
    d = ALPHA(args) * pow(2*M_PI*pow((T(args)-THETA(args)),3), -.5) * 
        exp( -pow(ALPHA(args)-NU(args)*(T(args)-THETA(args)),2) / (2*(T(args)-THETA(args))));
  
  return(d);
}


}

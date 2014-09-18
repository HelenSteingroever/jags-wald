#include <config.h>
#include "DWaldGamma.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <JRmath.h>
#include <cmath>

#include <gsl/gsl_specfunc.h>

using std::vector;

#define LAMBDA(par) (*par[0])
#define ALPHA(par) (*par[1])
#define TAU(par) (*par[2])
#define KAPPA(par) (*par[3])

namespace wald {

DWaldGamma::DWaldGamma() : ScalarDist("dwald_gamma", 4, DIST_PROPORTION)
{}

bool DWaldGamma::checkParameterValue (vector<double const *> const &parameters) const
{
    return  (true);
}

double DWaldGamma::r_gamma(double x) const
{
  return exp(gamma(x));
}

double DWaldGamma::LaguerreL(double n, double a, double x) const
{
  double Lf;

  Lf = 1/((n+a+1)*(r_gamma(1+a)*r_gamma(n+1)/r_gamma((1+a)+(n+1)))) * gsl_sf_hyperg_1F1(-n,a+1,x); 

  return Lf;
}

double DWaldGamma::erf(double x) const
{
  double e;
  
  e = 2 * pnorm(x * sqrt(2), 0,1, 1,0) - 1;
  
  return e;
}

double DWaldGamma::dwald_gamma(double t, vector<double const *> const &parameters) const
{
  double lambda = LAMBDA(parameters); 
  double alpha = ALPHA(parameters); 
  double tau = TAU(parameters); 
  double kappa = KAPPA(parameters);
  
  double d;
  
  if (t < 0.01 || tau >= 2 || kappa <= 0.2)
  {
    d = 0;  
  }
  else if (lambda == 1 && alpha == 1 && tau == 1 && kappa == 1 )
  {
    d = exp(- 1 / (2*t)) / (2 * pow(t,2));
  }
  else if (lambda == 1 && tau == 1) 
  {
    d = alpha*exp(-(2*alpha*kappa-1)/(2*t*pow(kappa,2)))* (1 + erf( (alpha*kappa-1) / (kappa*sqrt(2*t)) )) / (2*pow(t,2)*kappa);
  }
  else 
  {  
    d = -sqrt(0.2e1) * pow(0.2e1, tau / 0.2e1 + 0.1e1 / 0.2e1) * alpha * exp(-0.1e1 / t * alpha * alpha * lambda / 0.2e1) * pow(kappa, -tau - 0.3e1) * pow(lambda, -tau / 0.2e1 - 0.1e1) * pow(t, -tau / 0.2e1 - 0.3e1) * 0.3141592654e1 * (-sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.5e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * alpha * alpha * pow(kappa, 0.3e1) + sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.5e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * alpha * alpha * pow(kappa, 0.3e1) - sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.3e1 / 0.2e1) * pow(t, 0.3e1 / 0.2e1) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * pow(kappa, 0.3e1) + 0.2e1 * sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.3e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * alpha * kappa * kappa - 0.2e1 * sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(lambda, 0.3e1 / 0.2e1) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * alpha * kappa * kappa - sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * sqrt(lambda) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * kappa + sqrt(0.2e1) * LaguerreL(-tau / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * sqrt(lambda) * sqrt(t) * cos(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.2e1) * kappa + sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(alpha, 0.3e1) * pow(kappa, 0.3e1) * pow(lambda, 0.3e1) - sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * pow(alpha, 0.3e1) * pow(kappa, 0.3e1) * pow(lambda, 0.3e1) + sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * pow(kappa, 0.3e1) * lambda * lambda * tau * t - sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * pow(kappa, 0.3e1) * lambda * lambda * t - 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * alpha * kappa * kappa * lambda * lambda + 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * alpha * kappa * kappa * lambda * lambda - sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * kappa * kappa * lambda * tau * t + sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * kappa * kappa * lambda * t + 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * kappa * lambda - 0.3e1 * sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) * alpha * kappa * lambda - sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.1e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1) + sin(0.3141592654e1 * tau / 0.2e1) * r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) * LaguerreL(-tau / 0.2e1 + 0.1e1 / 0.2e1, 0.3e1 / 0.2e1, 0.1e1 / lambda * pow(alpha * kappa * lambda - 0.1e1, 0.2e1) * pow(kappa, -0.2e1) / t / 0.2e1)) / r_gamma(tau) / sin(0.3141592654e1 * tau / 0.2e1) / r_gamma(-tau / 0.2e1 + 0.3e1 / 0.2e1) / cos(0.3141592654e1 * tau / 0.2e1) / r_gamma(-tau / 0.2e1 + 0.2e1) / 0.16e2;
  }   

  if (d<0) d = 0;

  return d;
}

double DWaldGamma::logDensity(double x, PDFType type,
       vector<double const *> const &parameters,
       double const *lbound, double const *ubound) const 
{
    double d = 0;

    d = dwald_gamma(x , parameters);

    return d == 0 ? JAGS_NEGINF : log(d);
}

double DWaldGamma::randomSample(vector<double const *> const &parameters, 
         double const *lbound, double const *ubound,
         RNG *rng) const
{
    return fabs( rng->uniform() );
}

double DWaldGamma::typicalValue(vector<double const *> const &parameters,
         double const *lbound, double const *ubound) const
{
    return 0.5;
}

bool DWaldGamma::isDiscreteValued(vector<bool> const &mask) const
{
    return true;
}

}

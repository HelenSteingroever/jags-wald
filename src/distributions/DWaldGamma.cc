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
  // calculates log density
  double lambda = LAMBDA(parameters); 
  double alpha = ALPHA(parameters); 
  double tau = TAU(parameters); 
  double kappa = KAPPA(parameters);
  
  double d;
    
  if (alpha == 1 && tau == 1 && kappa == 1 )
  {
      d = (- 1. / (2.*t)) - log(2.) -2*log(t);
  }
  else if (tau == 1)
  {
      d = log(alpha) - ((2.*alpha*kappa-1.)/(2.*t*pow(kappa,2))) + 
      log(1. + erf( (alpha*kappa-1.) / (kappa*sqrt(2.*t)) )) -
      log(2.) - 2*log(t) - log(kappa);
  }
  else
  {
      double L1; 
      double L2; 
      double L3; 
      double L4; 
      double C1; 
      double C2; 
                                                 
      L1 = LaguerreL(-(.5)*tau+.5, .5, (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
      
      L2 = LaguerreL(-(.5)*tau+.5, 3./2., (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
      
      L3 = LaguerreL(-(.5)*tau, .5, (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
      
      L4 = LaguerreL(-(.5)*tau, 3./2., (.5)*pow(alpha*kappa-1.,2)/(pow(kappa,2)*t));
      
      C1 = sin((.5)*M_PI*tau)*r_gamma(-(.5)*tau+3./2.);
      
      C2 = sqrt(2.)*pow(alpha,3)*pow(kappa,3)*sqrt(t);
      
      d = log(-(1.0/16.0)*pow(2.,((.5)*tau+.5))*alpha*exp(-(.5)*pow(alpha,2)/t)*pow(kappa,(-tau-3.))*
      pow(t,(-(.5)*tau-7.0/2.0))*M_PI*
      
      (
       
       C1*
       L1*
       C2 +
       
       C1*
       L1*
       sqrt(2.)*alpha*pow(kappa,3)*tau*pow(t,(3./2.)) -
       
       C1*
       L2*
       C2 -
       
       C1*
       L1*
       sqrt(2.)*alpha*pow(kappa,3)*pow(t,(3./2.)) -
       
       3.*C1*
       L1*
       sqrt(2.)*pow(alpha,2)*pow(kappa,2)*sqrt(t) -
       
       C1*
       L1*
       sqrt(2.)*pow(kappa,2)*tau*pow(t,(3./2.)) +
       
       3.*C1*
       L2*
       sqrt(2.)*pow(alpha,2)*pow(kappa,2)*sqrt(t)-
       
       2.*
       L3*
       cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*pow(alpha,2)*pow(kappa,3)*t +
       
       2.*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*
       L4*
       pow(alpha,2)*pow(kappa,3)*t +
       
       C1*
       L1*
       sqrt(2.)*pow(kappa,2)*pow(t,(3./2.)) -
       
       2.*
       L3*
       cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*pow(kappa,3)*pow(t,2) +
       
       3.*C1*
       L1*
       sqrt(2.)*alpha*kappa*sqrt(t) -
       
       3.*C1*
       L2*
       sqrt(2.)*alpha*kappa*sqrt(t) +
       
       4.*
       L3*
       cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*alpha*pow(kappa,2)*t -
       
       4.*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*
       L4*
       alpha*pow(kappa,2)*t -
       
       C1*
       L1*
       sqrt(2.)*sqrt(t) +
       
       C1*
       L2*
       sqrt(2.)*sqrt(t) -
       
       2.*
       L3*
       cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*kappa*t +
       
       2.*cos((.5)*M_PI*tau)*r_gamma(-(.5)*tau+2.)*
       L4*
       kappa*t
       
       ))-
      
      log(r_gamma(tau)) -log(C1) - log(cos((.5)*M_PI*tau)) - log(r_gamma(-(.5)*tau+2.));
  }
  return(d);
}

double DWaldGamma::logDensity(double x, PDFType type,
       vector<double const *> const &parameters,
       double const *lbound, double const *ubound) const 
{
    double logd = 0;

    logd = dwald_gamma(x , parameters);

    return logd;
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

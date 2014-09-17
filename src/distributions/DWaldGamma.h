#ifndef DWALDGAMMA_H_
#define DWALDGAMMA_H_
#include <distribution/ScalarDist.h> // JAGS scalar distribution base class

namespace wald {

class DWaldGamma : public ScalarDist // scalar distribution class
{
  public:
    DWaldGamma(); // constructor
    double r_gamma(double x) const;
    double LaguerreL(double n, double a, double x) const;
    double erf(double x) const;
    double dwald_gamma(double t, std::vector<double const *> const &parameters) const;
    double logDensity(double x, PDFType type,
                      std::vector<double const *> const &parameters,
                      double const *lower, double const *upper) const;
    double randomSample(std::vector<double const *> const &parameters,
                        double const *lower, double const *upper,
                        RNG *rng) const;
    double typicalValue(std::vector<double const *> const &parameters,
                        double const *lower, double const *upper) const;
    /** Checks that p lies in the open interval (0,1) */
    bool checkParameterValue(std::vector<double const *> const &parameters) const;
    /** Distribution is discrete valued */
    bool isDiscreteValued(std::vector<bool> const &mask) const;
};

}
#endif /* DWALDGAMMA_H_ */

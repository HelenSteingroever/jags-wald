#ifndef DSWALDFUNC_H_
#define DSWALDFUNC_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace wald {

class DSWaldfunc : public ScalarFunction 
{
  public:
    DSWaldfunc();

    bool checkParameterValue(std::vector<double const *> const &args) const;
    double evaluate(std::vector<double const *> const &args) const;
};

} // namespace wald
} // namespace jags

#endif /* DSWALDFUNC_H_ */

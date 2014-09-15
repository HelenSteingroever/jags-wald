#ifndef BERN_FUNC_H_
#define BERN_FUNC_H_

#include <function/ScalarFunction.h>

namespace bernoulli {

class LogBernFun : public ScalarFunction 
{
  public:
    LogBernFun();

    bool checkParameterValue(std::vector<double const *> const &args) const;
    double evaluate(std::vector<double const *> const &args) const;
};

}

#endif /* BERN_FUNC_H_ */

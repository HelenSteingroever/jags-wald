#include <Module.h> // JAGS module base class
#include <distributions/DBern.h> // Bernoulli distribution class
#include <functions/LogBernFun.h> // Bernoulli log function class

namespace bernoulli { // module namespace

class BERNModule : public Module { // module class
  public:
    BERNModule(); // constructor
    ~BERNModule(); // destructor
};

BERNModule::BERNModule() : Module("bernoulli") {
  insert(new DBern); // inherited function to load JAGS objects
  insert(new LogBernFun); // inherited function to load JAGS objects
}
BERNModule::~BERNModule() {
  std::vector<Distribution*> const &dvec = distributions();
  for (unsigned int i = 0; i < dvec.size(); ++i) {
    delete dvec[i];
  } // deletes instantiated distribution objects
}

}

bernoulli::BERNModule _bernoulli_module;

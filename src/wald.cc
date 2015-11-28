#include <Module.h> // JAGS module base class
#include <distributions/DWald.h> // Wald (Inverse Gaussian) distribution class
#include <distributions/DSWald.h> // Shifted Wald distribution class
#include <distributions/PSWald.h> // Shifted Wald cumulative distribution class
#include <distributions/DWaldTrunc.h> // Driftrate as Truncated - Wald Normal distribution class
#include <distributions/DWaldGamma.h> // Driftrate as Gamma - Wald distribution class
#include <functions/DSWaldfunc.h> // 
#include <functions/PSWaldfunc.h> // 

namespace jags { // jags namespace
namespace wald { // module namespace

class WALDModule : public Module { // module class
  public:
    WALDModule(); // constructor
    ~WALDModule(); // destructor
};

WALDModule::WALDModule() : Module("wald") {
  // load distributions
  insert(new DWald); // inherited function to load JAGS objects
  insert(new DSWald); // inherited function to load JAGS objects
  insert(new PSWald); // inherited function to load JAGS objects
  insert(new DWaldTrunc); // inherited function to load JAGS objects
  insert(new DWaldGamma); // inherited function to load JAGS objects

  // load functions
  insert(new DSWaldfunc);
  insert(new PSWaldfunc);
}
WALDModule::~WALDModule() {
  std::vector<Distribution*> const &dvec = distributions();
  for (unsigned int i = 0; i < dvec.size(); ++i) {
    delete dvec[i];
  } // deletes instantiated distribution objects

  std::vector<Function*> const &fvec = functions();
  for (unsigned int i = 0; i < fvec.size(); ++i) {
    delete fvec[i];
  }
}

} // namespace wald
} // namespace jags

wald::WALDModule _wald_module;

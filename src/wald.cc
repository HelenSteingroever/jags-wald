#include <module/Module.h> // JAGS module base class
#include <distributions/DSWald.h> // Shifted Wald distribution class
#include <distributions/DSWTN.h> // Driftrate as Truncated - Wald Normal distribution class
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
  insert(new DSWTN); // inherited function to load JAGS objects
  insert(new DSWald); // inherited function to load JAGS objects

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

jags::wald::WALDModule _wald_module;

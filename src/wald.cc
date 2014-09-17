#include <Module.h> // JAGS module base class
#include <distributions/DWald.h> // Wald distribution class
#include <distributions/DWaldGamma.h> // Wald distribution class

namespace wald { // module namespace

class WALDModule : public Module { // module class
  public:
    WALDModule(); // constructor
    ~WALDModule(); // destructor
};

WALDModule::WALDModule() : Module("wald") {
  insert(new DWald); // inherited function to load JAGS objects
  insert(new DWaldGamma); // inherited function to load JAGS objects
}
WALDModule::~WALDModule() {
  std::vector<Distribution*> const &dvec = distributions();
  for (unsigned int i = 0; i < dvec.size(); ++i) {
    delete dvec[i];
  } // deletes instantiated distribution objects
}

}

wald::WALDModule _wald_module;

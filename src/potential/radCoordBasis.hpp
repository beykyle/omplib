#include "util/constants.hpp"

namespace omplib {

enum class Coord : bool {
  local,
  nonlocal
};

class LocalPotential {
public:
  constexpr static Coord coord = Coord::local;
  double operator () (double r) const;
};

class WoodSaxon : public LocalPotential {
};

class DerivWoodSaxon : public LocalPotential {
};

class Minnesota : public LocalPotential {
};

class Thompson : public LocalPotential {
};

class NonLocalPotential {
public:
  constexpr static Coord coord = Coord::nonlocal;
  double operator () (double r, double rp) const;
};

template<class LocalPotential>
class PeryBuck : NonLocalPotential {
  double operator () (double r, double rp) const;
private:
  LocalPotential potential;
};

}

//
// Scientific sanitizer checks on molecules (see RDGeneral/Sanitizer.h).
// Observation only; every function is a no-op unless RDKIT_SANITIZE_LOG is set
// and never modifies its arguments.
//
#include <RDGeneral/export.h>
#ifndef RD_SANITIZERMOL_H
#define RD_SANITIZERMOL_H

#include <RDGeneral/Sanitizer.h>
#include <RDGeneral/types.h>
#include <exception>
#include <map>
#include <string>
#include <vector>

namespace RDKit {
class ROMol;
namespace Sanitizer {

//! element counts (explicit atoms, plus implicit/explicit-count hydrogens
//! when `countImplicitH`) and net formal charge
struct Composition {
  bool valid = false;
  std::map<int, long> n;
  long charge = 0;
};
RDKIT_GRAPHMOL_EXPORT Composition compose(const ROMol &mol,
                                          bool countImplicitH);
//! RDK-MO-001/002/003, RDK-MS-001/002/004/005/006
RDKIT_GRAPHMOL_EXPORT void checkSameComposition(const char *id,
                                                const Composition &before,
                                                const ROMol &after,
                                                bool countImplicitH,
                                                const char *where);
//! RDK-MS-003
RDKIT_GRAPHMOL_EXPORT void checkProtonBookkeeping(const Composition &before,
                                                  const ROMol &after);

struct Valences {
  bool valid = false;
  std::vector<unsigned int> v;
};
//! RDK-MO-004
RDKIT_GRAPHMOL_EXPORT Valences snapValences(const ROMol &mol);
RDKIT_GRAPHMOL_EXPORT void checkKekulizeValences(const Valences &before,
                                                 const ROMol &after);

//! RDK-MO-005/006
RDKIT_GRAPHMOL_EXPORT void checkSSSR(const ROMol &mol,
                                     const VECT_INT_VECT &rings, int count,
                                     bool includeDative, bool includeHBonds);

//! RDK-SI-001
RDKIT_GRAPHMOL_EXPORT void checkSameConstitution(const ROMol &ref,
                                                 const ROMol &iso);

//! RAII: composition snapshot at construction, verdict at normal exit
class CompositionGuard {
 public:
  CompositionGuard(const char *id, const ROMol &mol, bool countImplicitH,
                   const char *where, bool active = true)
      : d_id(id), d_mol(mol), d_h(countImplicitH), d_where(where) {
    if (active && enabled()) {
      d_before = compose(mol, countImplicitH);
    }
  }
  ~CompositionGuard() {
    if (d_before.valid && std::uncaught_exceptions() == 0) {
      checkSameComposition(d_id, d_before, d_mol, d_h, d_where);
    }
  }
  CompositionGuard(const CompositionGuard &) = delete;
  CompositionGuard &operator=(const CompositionGuard &) = delete;

 private:
  const char *d_id;
  const ROMol &d_mol;
  bool d_h;
  const char *d_where;
  Composition d_before;
};

//! RAII: proton-bookkeeping snapshot (RDK-MS-003)
class ProtonGuard {
 public:
  explicit ProtonGuard(const ROMol &mol) : d_mol(mol) {
    if (enabled()) {
      d_before = compose(mol, true);
    }
  }
  ~ProtonGuard() {
    if (d_before.valid && std::uncaught_exceptions() == 0) {
      checkProtonBookkeeping(d_before, d_mol);
    }
  }
  ProtonGuard(const ProtonGuard &) = delete;
  ProtonGuard &operator=(const ProtonGuard &) = delete;

 private:
  const ROMol &d_mol;
  Composition d_before;
};

//! RAII: valence snapshot at construction, verdict at normal exit
class KekulizeGuard {
 public:
  explicit KekulizeGuard(const ROMol &mol) : d_mol(mol) {
    d_before = snapValences(mol);
  }
  ~KekulizeGuard() {
    if (d_before.valid && std::uncaught_exceptions() == 0) {
      checkKekulizeValences(d_before, d_mol);
    }
  }
  KekulizeGuard(const KekulizeGuard &) = delete;
  KekulizeGuard &operator=(const KekulizeGuard &) = delete;

 private:
  const ROMol &d_mol;
  Valences d_before;
};
}  // namespace Sanitizer
}  // namespace RDKit
#endif

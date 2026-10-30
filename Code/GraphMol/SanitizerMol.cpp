//
// Scientific sanitizer checks on molecules (see SanitizerMol.h).
//
#include "SanitizerMol.h"
#include <GraphMol/RDKitBase.h>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <sstream>

namespace RDKit {
namespace Sanitizer {
namespace {
std::string describe(const Composition &c) {
  std::ostringstream os;
  for (const auto &kv : c.n) {
    os << kv.first << ":" << kv.second << " ";
  }
  os << "q=" << c.charge;
  return os.str();
}
}  // namespace

Composition compose(const ROMol &mol, bool countImplicitH) {
  Composition c;
  try {
    for (const auto atom : mol.atoms()) {
      c.n[atom->getAtomicNum()] += 1;
      c.charge += atom->getFormalCharge();
      if (countImplicitH) {
        c.n[1] += atom->getNumExplicitHs() + atom->getNumImplicitHs();
      }
    }
  } catch (...) {
    return Composition();
  }
  for (auto it = c.n.begin(); it != c.n.end();) {
    it = it->second == 0 ? c.n.erase(it) : std::next(it);
  }
  c.valid = true;
  return c;
}

void checkSameComposition(const char *id, const Composition &before,
                          const ROMol &after, bool countImplicitH,
                          const char *where) {
  if (!enabled() || !before.valid) {
    return;
  }
  const Composition a = compose(after, countImplicitH);
  if (!a.valid) {
    return;
  }
  if (a.n != before.n || a.charge != before.charge) {
    report(id, std::string(where) + ": before [" + describe(before) +
                   "] after [" + describe(a) + "]");
  }
}

void checkProtonBookkeeping(const Composition &before, const ROMol &after) {
  if (!enabled() || !before.valid) {
    return;
  }
  const Composition a = compose(after, true);
  if (!a.valid) {
    return;
  }
  auto heavy = [](std::map<int, long> m) {
    m.erase(1);
    return m;
  };
  auto nH = [](const Composition &c) {
    auto it = c.n.find(1);
    return it == c.n.end() ? 0L : it->second;
  };
  if (heavy(a.n) != heavy(before.n) ||
      a.charge - nH(a) != before.charge - nH(before)) {
    report("RDK-MS-003", "before [" + describe(before) + "] after [" +
                             describe(a) + "]");
  }
}

Valences snapValences(const ROMol &mol) {
  Valences v;
  if (!enabled()) {
    return v;
  }
  try {
    for (const auto atom : mol.atoms()) {
      if (atom->needsUpdatePropertyCache()) {
        return Valences();
      }
      v.v.push_back(atom->getTotalValence());
    }
  } catch (...) {
    return Valences();
  }
  v.valid = true;
  return v;
}

void checkKekulizeValences(const Valences &before, const ROMol &after) {
  if (!enabled() || !before.valid || before.v.size() != after.getNumAtoms()) {
    return;
  }
  for (const auto atom : after.atoms()) {
    double sum = 0.0;
    for (const auto bond : after.atomBonds(atom)) {
      sum += bond->getValenceContrib(atom);
    }
    // total valence = bond valence plus every hydrogen not in the graph
    try {
      sum += atom->getNumExplicitHs() + atom->getNumImplicitHs();
    } catch (...) {
      return;
    }
    if (std::fabs(sum - before.v[atom->getIdx()]) > 1e-9) {
      report("RDK-MO-004", "atom " + std::to_string(atom->getIdx()) +
                               " valence before " +
                               std::to_string(before.v[atom->getIdx()]) +
                               " after " + std::to_string(sum));
      return;
    }
  }
}

void checkSSSR(const ROMol &mol, const VECT_INT_VECT &rings, int count,
               bool includeDative, bool includeHBonds) {
  if (!enabled()) {
    return;
  }
  const unsigned int n = mol.getNumAtoms();
  std::vector<unsigned int> parent(n);
  std::iota(parent.begin(), parent.end(), 0u);
  auto find = [&](unsigned int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  auto counted = [&](const Bond *b) {
    if (b->getBondType() == Bond::DATIVE ||
        b->getBondType() == Bond::DATIVEONE) {
      return includeDative;
    }
    if (b->getBondType() == Bond::HYDROGEN) {
      return includeHBonds;
    }
    return true;
  };
  long E = 0;
  for (const auto b : mol.bonds()) {
    if (!counted(b)) {
      continue;
    }
    ++E;
    parent[find(b->getBeginAtomIdx())] = find(b->getEndAtomIdx());
  }
  long C = 0;
  for (unsigned int i = 0; i < n; ++i) {
    C += find(i) == i;
  }
  const long cyclo = E - static_cast<long>(n) + C;
  if (count != cyclo) {
    report("RDK-MO-005", "SSSR count " + std::to_string(count) +
                             " cyclomatic number " + std::to_string(cyclo));
  }
  for (const auto &ring : rings) {
    bool ok = ring.size() >= 3;
    std::vector<int> sorted(ring);
    std::sort(sorted.begin(), sorted.end());
    ok = ok && std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
    for (size_t i = 0; ok && i < ring.size(); ++i) {
      const int a = ring[i], b = ring[(i + 1) % ring.size()];
      if (a < 0 || b < 0 || a >= static_cast<int>(n) ||
          b >= static_cast<int>(n)) {
        ok = false;
        break;
      }
      const Bond *bd = mol.getBondBetweenAtoms(a, b);
      ok = bd != nullptr && counted(bd);
    }
    if (!ok) {
      std::string s;
      for (int a : ring) {
        s += std::to_string(a) + " ";
      }
      report("RDK-MO-006", "ring is not a simple bonded cycle: " + s);
      return;
    }
  }
}

void checkSameConstitution(const ROMol &ref, const ROMol &iso) {
  if (!enabled()) {
    return;
  }
  bool same = ref.getNumAtoms() == iso.getNumAtoms() &&
              ref.getNumBonds() == iso.getNumBonds();
  for (unsigned int i = 0; same && i < ref.getNumAtoms(); ++i) {
    const auto a = ref.getAtomWithIdx(i), b = iso.getAtomWithIdx(i);
    same = a->getAtomicNum() == b->getAtomicNum() &&
           a->getFormalCharge() == b->getFormalCharge();
  }
  for (unsigned int i = 0; same && i < ref.getNumBonds(); ++i) {
    const auto a = ref.getBondWithIdx(i), b = iso.getBondWithIdx(i);
    same = a->getBeginAtomIdx() == b->getBeginAtomIdx() &&
           a->getEndAtomIdx() == b->getEndAtomIdx() &&
           a->getBondType() == b->getBondType();
  }
  if (!same) {
    report("RDK-SI-001", "stereoisomer differs in atoms, charges, bonds or bond orders");
    return;
  }
  const auto c0 = compose(ref, true), c1 = compose(iso, true);
  if (c0.valid && c1.valid && (c0.n != c1.n || c0.charge != c1.charge)) {
    report("RDK-SI-001", "stereoisomer composition differs");
  }
}
}  // namespace Sanitizer
}  // namespace RDKit

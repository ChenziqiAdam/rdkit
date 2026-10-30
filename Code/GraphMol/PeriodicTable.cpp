//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "PeriodicTable.h"
#include <GraphMol/SanitizerMol.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <atomic>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <sstream>
#include <locale>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif

namespace RDKit {

namespace {
std::atomic<PeriodicTable *> ds_rawInstance{nullptr};
}  // namespace

class std::unique_ptr<PeriodicTable> PeriodicTable::ds_instance = nullptr;

// RDK-PT-001..005: consistency of the tabulated element data
static void sanitizeTables(const std::vector<atomicData> &byanum) {
  if (!Sanitizer::enabled()) {
    return;
  }
  for (const auto &el : byanum) {
    const std::string sym = el.Symbol();
    double sum = 0.0, wmass = 0.0, amax = 0.0;
    for (const auto &kv : el.d_isotopeInfoMap) {
      const double m = kv.second.first, x = kv.second.second;
      sum += x;
      wmass += x * m / 100.0;
      amax = std::max(amax, x);
      // bound nuclide (naturally occurring isotope): binding energy between
      // 0 and 9 MeV per nucleon, so A (1 - 0.0099) <= m <= Z m(1H) + N m(n)
      const double A = kv.first, Z = el.AtomicNum();
      if (x > 0.0 && Z > 0 &&
          (m < A * (1.0 - 0.0099) ||
           m > Z * 1.00782503223 + (A - Z) * 1.00866491595 + 1e-6)) {
        Sanitizer::report("RDK-PT-003", sym + " isotope " +
                                            std::to_string(kv.first) +
                                            " mass " + std::to_string(m));
      }
    }
    if (sum > 0.0) {
      if (std::fabs(sum - 100.0) > 0.02) {
        Sanitizer::report("RDK-PT-001",
                          sym + " abundances sum " + std::to_string(sum));
      } else if (std::fabs(wmass - el.Mass()) > 2e-3 * el.Mass()) {
        Sanitizer::report("RDK-PT-002", sym + " weighted " +
                                            std::to_string(wmass) + " tabulated " +
                                            std::to_string(el.Mass()));
      }
      const int ci = el.MostCommonIsotope();
      auto it = el.d_isotopeInfoMap.find(ci);
      if (ci > 0 && (it == el.d_isotopeInfoMap.end() ||
                     it->second.second < amax ||
                     std::fabs(it->second.first - el.MostCommonIsotopeMass()) >
                         1e-6)) {
        Sanitizer::report("RDK-PT-005", sym + " common isotope " +
                                            std::to_string(ci));
      }
    }
    if (el.Rcov() > 0.0 && el.Rvdw() > 0.0 && el.Rvdw() <= el.Rcov()) {
      Sanitizer::report("RDK-PT-004",
                        sym + " Rvdw " + std::to_string(el.Rvdw()) + " Rcov " +
                            std::to_string(el.Rcov()));
    }
  }
}

PeriodicTable::PeriodicTable() {
  // it is assumed that the atomic atomData string constains atoms
  // in sequence and no atoms are missing in between
  byanum.clear();
  byname.clear();

  boost::char_separator<char> eolSep("\n");
  tokenizer tokens(periodicTableAtomData, eolSep);
  for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
       ++token) {
    if (*token != " ") {
      atomicData adata(*token);
      std::string enam = adata.Symbol();
      byname[enam] = adata.AtomicNum();
      // there are, for backwards compatibility reasons, some duplicate rows for
      // atomic numbers in the atomic_data data structure. It's ok to have
      // multiple symbols map to the same atomic number (above), but we need to
      // be sure that we only store one entry per atomic number.
      // Note that this only works because the first atom in the adata list is
      // the dummy atom (atomic number 0). This was #2784
      if (rdcast<size_t>(adata.AtomicNum()) == byanum.size()) {
        byanum.push_back(adata);
      }
    }
  }

  unsigned int lidx = 0;
  std::istringstream istr;
  istr.imbue(std::locale("C"));
  while (isotopesAtomData[lidx] != "" && isotopesAtomData[lidx] != "EOS") {
    tokenizer lines(isotopesAtomData[lidx++], eolSep);
    boost::char_separator<char> spaceSep(" \t");
    for (tokenizer::iterator line = lines.begin(); line != lines.end();
         ++line) {
      if (*line != " ") {
        tokenizer tokens(*line, spaceSep);
        tokenizer::iterator token = tokens.begin();
        int anum;
        istr.clear();
        istr.str(*token);
        istr >> anum;
        atomicData &adata = byanum[anum];
        ++token;
        if (token == tokens.end()) {
          continue;
        }
        ++token;
        if (token == tokens.end()) {
          continue;
        }
        unsigned int isotope;
        istr.clear();
        istr.str(*token);
        istr >> isotope;
        ++token;
        if (token == tokens.end()) {
          continue;
        }
        double mass;
        istr.clear();
        istr.str(*token);
        istr >> mass;
        ++token;
        if (token == tokens.end()) {
          continue;
        }
        double abundance;
        istr.clear();
        istr.str(*token);
        istr >> abundance;
        adata.d_isotopeInfoMap[isotope] = std::make_pair(mass, abundance);
      }
    }
  }
  sanitizeTables(byanum);
}

void PeriodicTable::initInstance() {
  ds_instance = std::unique_ptr<PeriodicTable>(new PeriodicTable());
  ds_rawInstance.store(ds_instance.get(), std::memory_order_release);
}

PeriodicTable *PeriodicTable::getTable() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  auto *raw = ds_rawInstance.load(std::memory_order_acquire);
  if (raw) {
    return raw;
  }
  static std::once_flag pt_init_once;
  std::call_once(pt_init_once, initInstance);
#else
  if (!ds_instance) {
    initInstance();
  }
#endif
  return ds_rawInstance.load(std::memory_order_acquire);
}

}  // namespace RDKit

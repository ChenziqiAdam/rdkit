//
//  Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BoundsMatrix.h"
#include "TriangleSmooth.h"
#include <RDGeneral/Sanitizer.h>
#include <algorithm>
#include <string>

namespace DistGeom {
namespace {
// RDK-DG-001..003 observation, run on a successfully smoothed matrix
void sanitizeBounds(const BoundsMatrix &b, double tol) {
  if (!RDKit::Sanitizer::enabled()) {
    return;
  }
  const unsigned int n = b.numRows();
  if (n > 150) {
    return;
  }
  double maxU = 0.0;
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = i + 1; j < n; ++j) {
      maxU = std::max(maxU, b.getUpperBound(i, j));
    }
  }
  // tol is relative (the smoothing code raises U to L within tol * L)
  const double slack = 8 * 2.220446049250313e-16 * maxU + tol * maxU;
  bool u = false, l = false, o = false;
  std::string du, dl, dord;
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = i + 1; j < n; ++j) {
      if (!o && b.getLowerBound(i, j) > b.getUpperBound(i, j) + slack) {
        o = true;
        dord = "pair " + std::to_string(i) + "," + std::to_string(j);
      }
      for (unsigned int k = 0; k < n; ++k) {
        if (k == i || k == j) {
          continue;
        }
        if (!u && b.getUpperBound(i, j) >
                      b.getUpperBound(i, k) + b.getUpperBound(k, j) + slack) {
          u = true;
          du = "triple " + std::to_string(i) + "," + std::to_string(j) + "," +
               std::to_string(k);
        }
        if (!l && b.getLowerBound(i, j) <
                      b.getLowerBound(i, k) - b.getUpperBound(k, j) - slack) {
          l = true;
          dl = "triple " + std::to_string(i) + "," + std::to_string(j) + "," +
               std::to_string(k);
        }
      }
    }
  }
  if (u) RDKit::Sanitizer::report("RDK-DG-001", du);
  if (l) RDKit::Sanitizer::report("RDK-DG-002", dl);
  if (o) RDKit::Sanitizer::report("RDK-DG-003", dord);
}
}  // namespace

bool triangleSmoothBounds(BoundsMatPtr boundsMat, double tol) {
  return triangleSmoothBounds(boundsMat.get(), tol);
}
bool triangleSmoothBounds(BoundsMatrix *boundsMat, double tol) {
  auto npt = boundsMat->numRows();
  for (auto k = 0u; k < npt; k++) {
    for (auto i = 0u; i < npt - 1; i++) {
      if (i == k) {
        continue;
      }
      auto ii = i;
      auto ik = k;
      if (ii > ik) {
        std::swap(ii, ik);
      }

      const auto Uik = boundsMat->getValUnchecked(ii, ik);  // upper bound
      const auto Lik = boundsMat->getValUnchecked(ik, ii);  // lower bound
      for (auto j = i + 1; j < npt; j++) {
        if (j == k) {
          continue;
        }
        auto jj = j;
        auto jk = k;
        if (jj > jk) {
          std::swap(jj, jk);
        }
        const auto Ukj = boundsMat->getValUnchecked(jj, jk);  // upper bound
        const auto sumUikUkj = Uik + Ukj;
        if (boundsMat->getValUnchecked(i, j) > sumUikUkj) {
          // adjust the upper bound
          boundsMat->setValUnchecked(i, j, sumUikUkj);
        }

        const auto diffLikUjk = Lik - Ukj;
        const auto diffLjkUik = boundsMat->getValUnchecked(jk, jj) - Uik;
        if (boundsMat->getValUnchecked(j, i) < diffLikUjk) {
          // adjust the lower bound
          boundsMat->setValUnchecked(j, i, diffLikUjk);
        } else if (boundsMat->getValUnchecked(j, i) < diffLjkUik) {
          // adjust the lower bound
          boundsMat->setValUnchecked(j, i, diffLjkUik);
        }
        const auto lBound = boundsMat->getValUnchecked(j, i);
        const auto uBound = boundsMat->getValUnchecked(i, j);
        if (tol > 0. && (lBound - uBound) / lBound > 0. &&
            (lBound - uBound) / lBound < tol) {
          // adjust the upper bound
          boundsMat->setValUnchecked(i, j, lBound);
        } else if (lBound - uBound > 0.) {
          return false;
        }
      }
    }
  }
  sanitizeBounds(*boundsMat, tol);
  return true;
}
}  // namespace DistGeom

//
// Scientific sanitizer runtime. Inactive unless the environment variable
// RDKIT_SANITIZE_LOG names a file; then every alarm appends one JSON line.
// Sanitizers only observe: they never change results, raise errors or alter
// state. Checks on plain numeric arrays live here; checks on molecules live in
// GraphMol/SanitizerMol.h.
//
#include <RDGeneral/export.h>
#ifndef RD_SANITIZER_H
#define RD_SANITIZER_H

#include <string>
#include <vector>

namespace RDKit {
namespace Sanitizer {

//! true when RDKIT_SANITIZE_LOG is set (evaluated once)
RDKIT_RDGENERAL_EXPORT bool enabled();

//! calls made by library code on its own intermediate molecules, which a user
//! never sees; observation of user-visible canonical forms is skipped inside
RDKIT_RDGENERAL_EXPORT void enterInternal();
RDKIT_RDGENERAL_EXPORT void leaveInternal();
RDKIT_RDGENERAL_EXPORT bool insideInternal();
class InternalScope {
 public:
  InternalScope() { enterInternal(); }
  ~InternalScope() { leaveInternal(); }
  InternalScope(const InternalScope &) = delete;
  InternalScope &operator=(const InternalScope &) = delete;
};

//! append one alarm record {"id","detail"}; no-op when disabled
RDKIT_RDGENERAL_EXPORT void report(const char *id, const std::string &detail);

//! RDK-FF-001/002. pos and gtot hold `npts` 3D points; fscale / tscale are
//! sums over every force-field contribution of |g| and |r - c||g|.
RDKIT_RDGENERAL_EXPORT void checkForceBalance(unsigned int npts,
                                              const double *pos,
                                              const double *gtot, double fscale,
                                              double tscale);

//! RDK-FF-003
RDKIT_RDGENERAL_EXPORT void checkEnergyDescent(double e0, double e1);

//! RDK-AL-001..003; points are mapped pairs (probe, reference), `trans` is a
//! row-major 4x4 matrix acting on column vectors; `msd` the reported mean
//! squared deviation; only called without weights.
RDKIT_RDGENERAL_EXPORT void checkAlignment(
    const std::vector<std::vector<double>> &prb,
    const std::vector<std::vector<double>> &ref, const double *trans,
    bool reflect, double msd);

//! RDK-MT-001/002; moments ascending, axes column-major 3x3 eigenvectors
RDKIT_RDGENERAL_EXPORT void checkInertia(const double *moments,
                                         const double *axes);

//! RDK-DESC-001
RDKIT_RDGENERAL_EXPORT void checkNPR(double pm1, double pm2, double pm3);

//! RDK-GAST-001
RDKIT_RDGENERAL_EXPORT void checkChargeConservation(double qStart, double qEnd,
                                                    unsigned int natoms,
                                                    unsigned int nIter,
                                                    double qmax);

//! RDK-MT-003/004. pos: positions of the atoms entering the sums, after
//! canonicalization; `origScale` = sum |r_orig||r_new|, `origNorm` = sum |r_orig|.
//! `centered`: the center argument was absent.
RDKIT_RDGENERAL_EXPORT void checkCanonicalFrame(
    const std::vector<std::vector<double>> &pos, double origScale,
    double origNorm, bool centered);
}  // namespace Sanitizer
}  // namespace RDKit
#endif

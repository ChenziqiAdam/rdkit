//
// Scientific sanitizer runtime (see Sanitizer.h). Observation only.
//
#include "Sanitizer.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <sstream>

namespace RDKit {
namespace Sanitizer {
namespace {
constexpr double EPS = 2.220446049250313e-16;
const char *logPath() {
  static const char *p = std::getenv("RDKIT_SANITIZE_LOG");
  return p;
}
std::mutex &logMutex() {
  static std::mutex m;
  return m;
}
std::string jsonEscape(const std::string &s) {
  std::string out;
  for (char c : s) {
    if (c == '"' || c == '\\') {
      out += '\\';
      out += c;
    } else if (c == '\n') {
      out += "\\n";
    } else if (static_cast<unsigned char>(c) >= 0x20) {
      out += c;
    }
  }
  return out;
}
using V3 = std::array<double, 3>;
V3 cross(const V3 &a, const V3 &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}
double norm(const V3 &a) {
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}
std::string num(double v) {
  std::ostringstream os;
  os.precision(10);
  os << v;
  return os.str();
}
}  // namespace

bool enabled() {
  const char *p = logPath();
  return p != nullptr && *p != '\0';
}

namespace {
thread_local int internalDepth = 0;
}
void enterInternal() { ++internalDepth; }
void leaveInternal() { --internalDepth; }
bool insideInternal() { return internalDepth > 0; }

void report(const char *id, const std::string &detail) {
  if (!enabled()) {
    return;
  }
  std::lock_guard<std::mutex> lock(logMutex());
  static std::map<std::string, int> counts;
  if (++counts[id] > 50) {
    return;
  }
  if (FILE *f = std::fopen(logPath(), "a")) {
    // the scorer sets RDKIT_SANITIZE_KEY per run; records without it are forgeries
    static const char *key = std::getenv("RDKIT_SANITIZE_KEY");
    std::fprintf(f, "{\"id\":\"%s\",\"key\":\"%s\",\"detail\":\"%s\"}\n", id,
                 key ? jsonEscape(key).c_str() : "", jsonEscape(detail).c_str());
    std::fclose(f);
  }
}

void checkForceBalance(unsigned int npts, const double *pos, const double *gtot,
                       double fscale, double tscale) {
  if (!enabled() || npts == 0) {
    return;
  }
  // rounding of each contribution is eps times its own magnitude; the factor
  // 1e-9 (~4e6 eps) leaves room for the conditioning of angle terms
  const double rel = 1e-9;
  V3 sum = {0, 0, 0}, c = {0, 0, 0};
  for (unsigned int i = 0; i < npts; ++i) {
    for (int d = 0; d < 3; ++d) {
      sum[d] += gtot[3 * i + d];
      c[d] += pos[3 * i + d] / npts;
    }
  }
  if (norm(sum) > rel * fscale) {
    report("RDK-FF-001", "net force " + num(norm(sum)) + " scale " +
                             num(fscale));
  }
  V3 tq = {0, 0, 0};
  for (unsigned int i = 0; i < npts; ++i) {
    V3 r = {pos[3 * i] - c[0], pos[3 * i + 1] - c[1], pos[3 * i + 2] - c[2]};
    V3 g = {gtot[3 * i], gtot[3 * i + 1], gtot[3 * i + 2]};
    V3 t = cross(r, g);
    for (int d = 0; d < 3; ++d) {
      tq[d] += t[d];
    }
  }
  if (norm(tq) > rel * tscale) {
    report("RDK-FF-002", "net torque " + num(norm(tq)) + " scale " +
                             num(tscale));
  }
}

void checkEnergyDescent(double e0, double e1) {
  if (!enabled() || !std::isfinite(e0) || !std::isfinite(e1)) {
    return;
  }
  if (e1 > e0 + 1e3 * EPS * std::max(1.0, std::fabs(e0))) {
    report("RDK-FF-003", "E0 " + num(e0) + " E1 " + num(e1));
  }
}

void checkAlignment(const std::vector<std::vector<double>> &prb,
                    const std::vector<std::vector<double>> &ref,
                    const double *t, bool reflect, double msd) {
  if (!enabled() || prb.size() != ref.size() || prb.empty() ||
      !std::isfinite(msd)) {
    return;
  }
  // 3x3 block: orthogonal, determinant sign tied to `reflect`
  double maxdev = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double s = 0.0;
      for (int k = 0; k < 3; ++k) {
        s += t[4 * k + i] * t[4 * k + j];
      }
      maxdev = std::max(maxdev, std::fabs(s - (i == j ? 1.0 : 0.0)));
    }
  }
  double det = t[0] * (t[5] * t[10] - t[6] * t[9]) -
               t[1] * (t[4] * t[10] - t[6] * t[8]) +
               t[2] * (t[4] * t[9] - t[5] * t[8]);
  if (maxdev > 1e-9 || (!reflect && det < 0.0) || std::fabs(std::fabs(det) - 1.0) > 1e-9) {
    report("RDK-AL-001", "orthogonality deviation " + num(maxdev) + " det " +
                             num(det) + " reflect " + (reflect ? "1" : "0"));
  }
  const double n = static_cast<double>(prb.size());
  V3 cp = {0, 0, 0}, cr = {0, 0, 0};
  for (size_t i = 0; i < prb.size(); ++i) {
    for (int d = 0; d < 3; ++d) {
      cp[d] += prb[i][d] / n;
      cr[d] += ref[i][d] / n;
    }
  }
  double mp = 0.0, mr = 0.0;
  for (size_t i = 0; i < prb.size(); ++i) {
    for (int d = 0; d < 3; ++d) {
      mp += (prb[i][d] - cp[d]) * (prb[i][d] - cp[d]) / n;
      mr += (ref[i][d] - cr[d]) * (ref[i][d] - cr[d]) / n;
    }
  }
  // offset of the two frames also enters the cancellation of the eigenvalue
  double off = 0.0;
  for (int d = 0; d < 3; ++d) {
    off += (cp[d] - cr[d]) * (cp[d] - cr[d]);
  }
  const double tol = 256 * EPS * (mp + mr + off + 0.0);
  auto apply = [&](const V3 &p) {
    V3 q;
    for (int i = 0; i < 3; ++i) {
      q[i] = t[4 * i] * p[0] + t[4 * i + 1] * p[1] + t[4 * i + 2] * p[2] +
             t[4 * i + 3];
    }
    return q;
  };
  double msd2 = 0.0;
  for (size_t i = 0; i < prb.size(); ++i) {
    V3 q = apply({prb[i][0], prb[i][1], prb[i][2]});
    for (int d = 0; d < 3; ++d) {
      msd2 += (q[d] - ref[i][d]) * (q[d] - ref[i][d]) / n;
    }
  }
  if (std::fabs(msd - msd2) > tol) {
    report("RDK-AL-002", "reported msd " + num(msd) + " recomputed " + num(msd2));
  }
  if (reflect) {
    return;
  }
  // local optimality over 64 fixed-seed small rotations about the centroid
  // of the transformed probe
  std::vector<V3> q(prb.size());
  V3 cq = {0, 0, 0};
  for (size_t i = 0; i < prb.size(); ++i) {
    q[i] = apply({prb[i][0], prb[i][1], prb[i][2]});
    for (int d = 0; d < 3; ++d) {
      cq[d] += q[i][d] / n;
    }
  }
  unsigned long long s = 88172645463325252ULL;
  auto rnd = [&]() {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return (s >> 11) * (1.0 / 9007199254740992.0) * 2.0 - 1.0;
  };
  for (int trial = 0; trial < 64; ++trial) {
    V3 ax = {rnd(), rnd(), rnd()};
    double an = norm(ax);
    if (an < 1e-3) {
      continue;
    }
    for (auto &x : ax) {
      x /= an;
    }
    double ang = 0.5 * rnd();
    double ca = std::cos(ang), sa = std::sin(ang);
    double m2 = 0.0;
    for (size_t i = 0; i < q.size(); ++i) {
      V3 p = {q[i][0] - cq[0], q[i][1] - cq[1], q[i][2] - cq[2]};
      V3 cr3 = cross(ax, p);
      double dot = ax[0] * p[0] + ax[1] * p[1] + ax[2] * p[2];
      for (int d = 0; d < 3; ++d) {
        double rot = p[d] * ca + cr3[d] * sa + ax[d] * dot * (1 - ca) + cq[d];
        m2 += (rot - ref[i][d]) * (rot - ref[i][d]) / n;
      }
    }
    if (m2 < msd2 - tol) {
      report("RDK-AL-003", "rotated msd " + num(m2) + " < aligned " + num(msd2));
      break;
    }
  }
}

void checkInertia(const double *m, const double *a) {
  if (!enabled()) {
    return;
  }
  double mx = std::max({std::fabs(m[0]), std::fabs(m[1]), std::fabs(m[2])});
  if (!std::isfinite(mx)) {
    return;
  }
  if (mx > 0.0) {
    double tol = 64 * EPS * mx;
    bool bad = m[0] < -tol || m[1] < -tol || m[2] < -tol;
    bad = bad || m[0] + m[1] - m[2] < -tol || m[0] + m[2] - m[1] < -tol ||
          m[1] + m[2] - m[0] < -tol;
    if (bad) {
      report("RDK-MT-001", "moments " + num(m[0]) + " " + num(m[1]) + " " +
                               num(m[2]));
    }
  }
  double dev = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double s = 0.0;
      for (int k = 0; k < 3; ++k) {
        s += a[3 * i + k] * a[3 * j + k];  // column-major: column i = a[3i..]
      }
      dev = std::max(dev, std::fabs(s - (i == j ? 1.0 : 0.0)));
    }
  }
  if (dev > 1e-12) {
    report("RDK-MT-002", "axes orthogonality deviation " + num(dev));
  }
}

void checkNPR(double p1, double p2, double p3) {
  if (!enabled() || !std::isfinite(p1 + p2 + p3) || p3 < 1e-8) {
    return;
  }
  double n1 = p1 / p3, n2 = p2 / p3, tol = 64 * EPS;
  if (n1 > n2 + tol || n2 > 1.0 + tol || n1 + n2 < 1.0 - tol) {
    report("RDK-DESC-001", "NPR1 " + num(n1) + " NPR2 " + num(n2));
  }
}

void checkChargeConservation(double q0, double q1, unsigned int natoms,
                             unsigned int nIter, double qmax) {
  if (!enabled() || !std::isfinite(q1)) {
    return;
  }
  double tol = 64 * EPS * (natoms + 1.0) * (nIter + 1.0) * std::max(1.0, qmax);
  if (std::fabs(q1 - q0) > tol) {
    report("RDK-GAST-001", "total charge start " + num(q0) + " end " + num(q1));
  }
}

void checkCanonicalFrame(const std::vector<std::vector<double>> &pos,
                         double origScale, double origNorm, bool centered) {
  if (!enabled() || pos.size() < 2) {
    return;
  }
  double m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  V3 s = {0, 0, 0};
  for (const auto &p : pos) {
    for (int i = 0; i < 3; ++i) {
      s[i] += p[i];
      for (int j = 0; j < 3; ++j) {
        m[i][j] += p[i] * p[j];
      }
    }
  }
  double tol = 1e3 * EPS * std::max(origScale, 1e-300);
  double off = std::max({std::fabs(m[0][1]), std::fabs(m[0][2]), std::fabs(m[1][2])});
  bool bad = off > tol || m[0][0] < m[1][1] - tol || m[1][1] < m[2][2] - tol;
  if (bad) {
    report("RDK-MT-003", "offdiag " + num(off) + " diag " + num(m[0][0]) + " " +
                             num(m[1][1]) + " " + num(m[2][2]));
  }
  if (centered && norm(s) > 1e3 * EPS * std::max(origNorm, 1e-300)) {
    report("RDK-MT-004", "centroid sum " + num(norm(s)));
  }
}
}  // namespace Sanitizer
}  // namespace RDKit

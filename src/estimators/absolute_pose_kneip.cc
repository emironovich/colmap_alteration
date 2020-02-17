// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "estimators/absolute_pose_kneip.h"
#include "estimators/absolute_pose.h"

#include "base/polynomial.h"
#include "estimators/utils.h"
#include "util/logging.h"
#include <iostream>

namespace colmap {
namespace {

Eigen::Vector3d LiftImagePoint(const Eigen::Vector2d& point) {
  return point.homogeneous() / std::sqrt(point.squaredNorm() + 1);
}

double sgn(double val) { return (double(0) < val) - (val < double(0)); }

}  // namespace

std::vector<P3PEstimatorKneip::M_t> P3PEstimatorKneip::Estimate(
    const std::vector<X_t>& points2D, const std::vector<Y_t>& points3D) {
  CHECK_EQ(points2D.size(), 3);
  CHECK_EQ(points3D.size(), 3);

  Eigen::Matrix3d f, P;
  P.col(0) = points3D[0];
  P.col(1) = points3D[1];
  P.col(2) = points3D[2];
  f.col(0) = LiftImagePoint(points2D[0]);
  f.col(1) = LiftImagePoint(points2D[1]);
  f.col(2) = LiftImagePoint(points2D[2]);

  //  P.col(0) << 0.547216, 0.138624, 0.149294;
  //  P.col(1) << 0.257508, 0.840717, 0.254282;
  //  P.col(2) << 0.814285, 0.243525, 0.929264;
  //
  //  f.col(0) << 0.508279, -0.628451, -0.588813;
  //  f.col(1) << -0.007459, 0.770069, -0.637917;
  //  f.col(2) << 0.905568, -0.063135, 0.419476;

  std::cout << "f:\n" << f << "\n";
  std::cout << "P:\n" << P << "\n";

  // Compute transformation matrix double and vector f_3^double.
  Eigen::Vector3d tx, ty, tz, f3_T;
  tx = f.col(0);
  tz = tx.cross(f.col(1));
  tz.normalize();
  ty = tz.cross(tx);
  Eigen::Matrix3d T;
  T << tx.transpose(), ty.transpose(), tz.transpose();
  f3_T = T * f.col(2);
  //  std::cout << "T:\n" << T << "\n";

  // Compute transformation matrix N and world point P_3^n.
  // TODO check nu exists by checking P1P2xP1P3 ~= 0
  Eigen::Vector3d nx, ny, nz, P3_n;
  nx = (P.col(1) - P.col(0));
  nx.normalize();
  nz = nx.cross(P.col(2) - P.col(0));
  nz.normalize();
  ny = nz.cross(nx);
  Eigen::Matrix3d N;
  N << nx.transpose(), ny.transpose(), nz.transpose();
  P3_n = N * (P.col(2) - P.col(0));
  //  std::cout << "N:\n" << N << "\n";
  //  std::cout << "P3_n:\n" << P3_n << "\n";

  // Extract p1 and p2 from P_3^n.
  // P3_n(3) = 0 check??????
  double p1 = P3_n(0);
  double p2 = P3_n(1);

  // Compute d12 and b.
  double d12 = (P.col(1) - P.col(0)).norm();
  double cos_beta = (f.col(0)).dot(f.col(1));
  double b = sgn(cos_beta) * (std::sqrt(1 / (1 - std::pow(cos_beta, 2)) - 1));

  //  std::cout << "d12:\n" << d12 << "\n";
  //  std::cout << "b:\n" << b << "\n";

  // Compute phi1 and phi2.
  double phi1 = f3_T(0) / f3_T(2);
  double phi2 = f3_T(1) / f3_T(2);

  // Compute factors a0 .... a4 of the polynomial.
  Eigen::Matrix<double, 1, 5> A;
  A(0) = -std::pow(phi2, 2) * std::pow(p2, 4) -
         std::pow(phi1, 2) * std::pow(p2, 4) - std::pow(p2, 4);
  A(1) = 2 * std::pow(p2, 3) * d12 * b +
         2 * std::pow(phi2, 2) * std::pow(p2, 3) * d12 * b -
         2 * phi1 * phi2 * std::pow(p2, 3) * d12;
  A(2) =
      -std::pow(phi2, 2) * std::pow(p1, 2) * std::pow(p2, 2) -
      std::pow(phi2, 2) * std::pow(p2, 2) * std::pow(d12, 2) * std::pow(b, 2) -
      std::pow(phi2, 2) * std::pow(p2, 2) * std::pow(d12, 2) +
      std::pow(phi2, 2) * std::pow(p2, 4) +
      std::pow(phi1, 2) * std::pow(p2, 4) + 2 * p1 * std::pow(p2, 2) * d12 +
      2 * phi1 * phi2 * p1 * std::pow(p2, 2) * d12 * b -
      std::pow(phi1, 2) * std::pow(p1, 2) * std::pow(p2, 2) +
      2 * std::pow(phi2, 2) * p1 * std::pow(p2, 2) * d12 -
      std::pow(p2, 2) * std::pow(d12, 2) * std::pow(b, 2) -
      2 * std::pow(p1, 2) * std::pow(p2, 2);
  A(3) = 2 * std::pow(p1, 2) * p2 * d12 * b +
         2 * phi1 * phi2 * std::pow(p2, 3) * d12 -
         2 * std::pow(phi2, 2) * std::pow(p2, 3) * d12 * b -
         2 * p1 * p2 * std::pow(d12, 2) * b;
  A(4) =
      -2 * phi1 * phi2 * p1 * std::pow(p2, 2) * d12 * b +
      std::pow(phi2, 2) * std::pow(p2, 2) * std::pow(d12, 2) +
      2 * std::pow(p1, 3) * d12 - std::pow(p1, 2) * std::pow(d12, 2) +
      std::pow(phi2, 2) * std::pow(p1, 2) * std::pow(p2, 2) - std::pow(p1, 4) -
      2 * std::pow(phi2, 2) * p1 * std::pow(p2, 2) * d12 +
      std::pow(phi1, 2) * std::pow(p1, 2) * std::pow(p2, 2) +
      std::pow(phi2, 2) * std::pow(p2, 2) * std::pow(d12, 2) * std::pow(b, 2);

  //  std::cout << "A:\n" << A << "\n";
  // Find roots for cos(theta)
  Eigen::Matrix<std::complex<double>, 4, 1> cos_theta_complex;
  double p =
      (8 * A(0) * A(2) - 3 * std::pow(A(1), 2)) / (8 * std::pow(A(0), 2));
  double q = (std::pow(A(1), 3) - 4 * A(0) * A(1) * A(2) +
              8 * std::pow(A(0), 2) * A(3)) /
             (8 * std::pow(A(0), 3));

  double d0 = std::pow(A(2), 2) - 3 * A(1) * A(3) + 12 * A(0) * A(4);
  double d1 = 2 * std::pow(A(2), 3) - 9 * A(1) * A(2) * A(3) +
              27 * std::pow(A(1), 2) * A(4) + 27 * A(0) * std::pow(A(3), 2) -
              72 * A(0) * A(2) * A(4);

  std::complex<double> Q =
      std::cbrt((d1 + std::sqrt(std::pow(d1, 2) - 4 * std::pow(d0, 3))) / 2);
  std::complex<double> S =
      std::sqrt(-2 * p / 3 + 1 / (3 * A(0)) * (Q + d0 / Q)) / 2.0;

  cos_theta_complex(0) = -A(1) / (4 * A(0)) - S +
                         std::sqrt(-4.0 * std::pow(S, 2) - 2 * p + q / S) / 2.0;
  cos_theta_complex(1) = -A(1) / (4 * A(0)) - S -
                         std::sqrt(-4.0 * std::pow(S, 2) - 2 * p + q / S) / 2.0;
  cos_theta_complex(2) = -A(1) / (4 * A(0)) + S +
                         std::sqrt(-4.0 * std::pow(S, 2) - 2 * p - q / S) / 2.0;
  cos_theta_complex(3) = -A(1) / (4 * A(0)) + S -
                         std::sqrt(-4.0 * std::pow(S, 2) - 2 * p - q / S) / 2.0;

  Eigen::Matrix<std::complex<double>, 5, 4> mons;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 5; ++j) {
      mons(j, i) = std::pow(cos_theta_complex(i), 4 - j);
    }
  }
  // std::cout << A*mons << std::endl;
  double e = 1e-4;  // todo fix
  std::vector<M_t> res;
  for (int i = 0; i < 4; ++i) {
    if ((cos_theta_complex(i)).imag() < e) {
      double cos_theta = (cos_theta_complex(i)).real();
      //      std::cout << "cos_theta:\n" << cos_theta << "\n";
      // Find cot(alpha) for each cos(theta)
      double cot_alpha = (phi1 * p1 / phi2 + cos_theta * p2 - d12 * b) /
                         (phi1 * cos_theta * p2 / phi2 - p1 + d12);

      // Compute all necessary trigonometric forms of alpha and theta
      double sin_alpha = std::sqrt(1 / (1 + std::pow(cot_alpha, 2)));
      double cos_alpha =
          sgn(cot_alpha) * std::sqrt(1 / (1 + 1 / (std::pow(cot_alpha, 2))));
      double sin_theta = -sgn(f3_T(2)) * std::sqrt(1 - std::pow(cos_theta, 2));

      // Compute C^n and Q
      Eigen::Vector3d C_n;
      C_n << d12 * cos_alpha * (sin_alpha * b + cos_alpha),
          d12 * sin_alpha * cos_theta * (sin_alpha * b + cos_alpha),
          d12 * sin_alpha * sin_theta * (sin_alpha * b + cos_alpha);
      Eigen::Matrix3d Q_m;
      Q_m << -cos_alpha, -sin_alpha * cos_theta, -sin_alpha * sin_theta,
          sin_alpha, -cos_alpha * cos_theta, -cos_alpha * sin_theta, 0,
          -sin_theta, cos_theta;
      //
      //      std::cout << "C_n:\n" << C_n << "\n";
      //      std::cout << "Q:\n" << Q_m << "\n";
      // Compute C and R
      M_t ans;
      Eigen::Vector3d C = P.col(0) + N.transpose() * C_n;
      Eigen::Matrix3d R = N.transpose() * Q_m.transpose() * T;
      std::cout << "C:\n" << C << "\n";
      std::cout << "R:\n" << R << "\n";
      ans << R, -R * C;
      res.push_back(ans);
    }
  }
  return res;
}

}  // namespace colmap

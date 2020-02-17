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

#define TEST_NAME "base/absolute_pose"
#include "util/testing.h"

#include <Eigen/Core>

#include "base/pose.h"
#include "base/similarity_transform.h"
#include "estimators/absolute_pose.h"
//#include "estimators/absolute_pose_kneip.h" //todo return tests
#include "estimators/essential_matrix.h"
#include "optim/ransac.h"
#include "util/random.h"

using namespace colmap;

BOOST_AUTO_TEST_CASE(TestP3P) {
  SetPRNGSeed(0);

  std::vector<Eigen::Vector3d> points3D;
  points3D.emplace_back(1, 1, 1);
  points3D.emplace_back(0, 1, 1);
  points3D.emplace_back(3, 1.0, 4);
  points3D.emplace_back(3, 1.1, 4);
  points3D.emplace_back(3, 1.2, 4);
  points3D.emplace_back(3, 1.3, 4);
  points3D.emplace_back(3, 1.4, 4);
  points3D.emplace_back(2, 1, 7);

  auto points3D_faulty = points3D;
  for (size_t i = 0; i < points3D.size(); ++i) {
    points3D_faulty[i](0) = 20;
  }

  for (double qx = 0; qx < 1; qx += 0.2) {
    for (double tx = 0; tx < 1; tx += 0.1) {
      const SimilarityTransform3 orig_tform(1, Eigen::Vector4d(1, qx, 0, 0),
                                            Eigen::Vector3d(tx, 0, 0));

      // Project points to camera coordinate system.
      std::vector<Eigen::Vector2d> points2D;
      for (size_t i = 0; i < points3D.size(); ++i) {
        Eigen::Vector3d point3D_camera = points3D[i];
        orig_tform.TransformPoint(&point3D_camera);
        points2D.push_back(point3D_camera.hnormalized());
      }

      RANSACOptions options;
      options.max_error = 1e-5;
      RANSAC<P3PEstimator> ransac(options);
      const auto report = ransac.Estimate(points2D, points3D);

      BOOST_CHECK_EQUAL(report.success, true);

      // Test if correct transformation has been determined.
      const double matrix_diff =
          (orig_tform.Matrix().topLeftCorner<3, 4>() - report.model).norm();
      BOOST_CHECK(matrix_diff < 1e-2);

      // Test residuals of exact points.
      std::vector<double> residuals;
      P3PEstimator::Residuals(points2D, points3D, report.model, &residuals);
      for (size_t i = 0; i < residuals.size(); ++i) {
        BOOST_CHECK(residuals[i] < 1e-3);
      }

      // Test residuals of faulty points.
      P3PEstimator::Residuals(points2D, points3D_faulty, report.model,
                              &residuals);
      for (size_t i = 0; i < residuals.size(); ++i) {
        BOOST_CHECK(residuals[i] > 0.1);
      }
    }
  }
}

// BOOST_AUTO_TEST_CASE(TestP3PKneip) {
//  SetPRNGSeed(0);
//
//  std::vector<Eigen::Vector3d> points3D;
//  points3D.emplace_back(1, 1, 1);
//  points3D.emplace_back(0, 1, 1);
//  points3D.emplace_back(3, 1.0, 4);
//  points3D.emplace_back(3, 1.1, 4);
//  points3D.emplace_back(3, 1.2, 4);
//  points3D.emplace_back(3, 1.3, 4);
//  points3D.emplace_back(3, 1.4, 4);
//  points3D.emplace_back(2, 1, 7);
//
//  auto points3D_faulty = points3D;
//  for (size_t i = 0; i < points3D.size(); ++i) {
//    points3D_faulty[i](0) = 20;
//  }
//
//  for (double qx = 0; qx < 1; qx += 0.2) {
//    for (double tx = 0; tx < 1; tx += 0.1) {
//      const SimilarityTransform3 orig_tform(1, Eigen::Vector4d(1, qx, 0, 0),
//                                            Eigen::Vector3d(tx, 0, 0));
//
//      // Project points to camera coordinate system.
//      std::vector<Eigen::Vector2d> points2D;
//      for (size_t i = 0; i < points3D.size(); ++i) {
//        Eigen::Vector3d point3D_camera = points3D[i];
//        orig_tform.TransformPoint(&point3D_camera);
//        points2D.push_back(point3D_camera.hnormalized());
//      }
//
//      RANSACOptions options;
//      options.max_error = 1e-5;
//      RANSAC<P3PEstimatorKneip> ransac(options);
//      const auto report = ransac.Estimate(points2D, points3D);
//
//      BOOST_CHECK_EQUAL(report.success, true);
//
//      // Test if correct transformation has been determined.
//      const double matrix_diff =
//          (orig_tform.Matrix().topLeftCorner<3, 4>() - report.model).norm();
//      BOOST_CHECK(matrix_diff < 1e-2);
//
//      // Test residuals of exact points.
//      std::vector<double> residuals;
//      P3PEstimatorKneip::Residuals(points2D, points3D, report.model,
//                                   &residuals);
//      for (size_t i = 0; i < residuals.size(); ++i) {
//        BOOST_CHECK(residuals[i] < 1e-3);
//      }
//
//      // Test residuals of faulty points.
//      P3PEstimatorKneip::Residuals(points2D, points3D_faulty, report.model,
//                                   &residuals);
//      for (size_t i = 0; i < residuals.size(); ++i) {
//        BOOST_CHECK(residuals[i] > 0.1);
//      }
//    }
//  }
//}

// BOOST_AUTO_TEST_CASE(TestP3PKneip) {
//  int num_samples = 3;
//
//  std::random_device dev;
//  std::mt19937_64 generator(dev());  // Mersenne Twister generator
//  std::uniform_real_distribution<double> uniformDistribution(-1., 1.);
//  auto uniform = [&]() { return uniformDistribution(generator); };
//
//  // rotation
//  Eigen::Vector4d rQuat = Eigen::Vector4d::NullaryExpr(4, 1, uniform);
//  Eigen::Matrix3d R = QuaternionToRotationMatrix(rQuat);
//
//  // camera center
//  Eigen::Vector3d C = (Eigen::Vector3d::NullaryExpr(3, 1, uniform)) / 2;
//
//  std::vector<Eigen::Vector3d> P;
//  std::vector<Eigen::Vector2d> points2D;
//  // std::vector<Eigen::Vector3d> f;
//  //  for (size_t i = 0; i < num_samples; ++i) {
//  //    Eigen::Vector3d pnt;
//  //    P.push_back(Eigen::Vector3d::NullaryExpr(3, 1, uniform));
//  //    pnt = R.transpose()*(P.back() - C);
//  //    points2D.emplace_back(pnt.hnormalized());
//  //  }
//
//  Eigen::Matrix3d PP, ff;
//  PP.col(0) << 0.547216, 0.138624, 0.149294;
//  PP.col(1) << 0.257508, 0.840717, 0.254282;
//  PP.col(2) << 0.814285, 0.243525, 0.929264;
//
//  ff.col(0) << 0.508279, -0.628451, -0.588813;
//  ff.col(1) << -0.007459, 0.770069, -0.637917;
//  ff.col(2) << 0.905568, -0.063135, 0.419476;
//
//  for (size_t i = 0; i < num_samples; ++i) {
//    Eigen::Vector3d pnt;
//    P.push_back(PP.col(i));
//    pnt = ff.col(i);
//    points2D.emplace_back(pnt.hnormalized());
//  }
//
//  std::vector<P3PEstimatorKneip::M_t> res =
//      P3PEstimatorKneip::Estimate(points2D, P);
//
//  for (auto M : res) {
//    std::cout << M << std::endl;
//  }
//
//  double min_diff_C = 1e+7;
//  double min_diff_R = 1e+7;
//
//  Eigen::Matrix3d R_es;
//  Eigen::Vector3d T_es;
//
//  for (auto M : res) {
//    R_es = M.leftCols<3>();
//    T_es = M.rightCols<1>();
//    if ((-R * C - T_es).norm() < 1e-4) {
//      min_diff_C = (-R * C - T_es).norm();
//      min_diff_R = (R - R_es).norm();
//    }
//  }
//
//  BOOST_CHECK(min_diff_C < 1e-4);
//  BOOST_CHECK(min_diff_R < 1e-4);
//}

BOOST_AUTO_TEST_CASE(TestEPNP) {
  SetPRNGSeed(0);

  std::vector<Eigen::Vector3d> points3D;
  points3D.emplace_back(1, 1, 1);
  points3D.emplace_back(0, 1, 1);
  points3D.emplace_back(3, 1.0, 4);
  points3D.emplace_back(3, 1.1, 4);
  points3D.emplace_back(3, 1.2, 4);
  points3D.emplace_back(3, 1.3, 4);
  points3D.emplace_back(3, 1.4, 4);
  points3D.emplace_back(2, 1, 7);

  auto points3D_faulty = points3D;
  for (size_t i = 0; i < points3D.size(); ++i) {
    points3D_faulty[i](0) = 20;
  }

  for (double qx = 0; qx < 1; qx += 0.2) {
    for (double tx = 0; tx < 1; tx += 0.1) {
      const SimilarityTransform3 orig_tform(1, Eigen::Vector4d(1, qx, 0, 0),
                                            Eigen::Vector3d(tx, 0, 0));

      // Project points to camera coordinate system.
      std::vector<Eigen::Vector2d> points2D;
      for (size_t i = 0; i < points3D.size(); ++i) {
        Eigen::Vector3d point3D_camera = points3D[i];
        orig_tform.TransformPoint(&point3D_camera);
        points2D.push_back(point3D_camera.hnormalized());
      }

      RANSACOptions options;
      options.max_error = 1e-5;
      RANSAC<EPNPEstimator> ransac(options);
      const auto report = ransac.Estimate(points2D, points3D);

      BOOST_CHECK_EQUAL(report.success, true);

      // Test if correct transformation has been determined.
      const double matrix_diff =
          (orig_tform.Matrix().topLeftCorner<3, 4>() - report.model).norm();
      BOOST_CHECK(matrix_diff < 1e-3);

      // Test residuals of exact points.
      std::vector<double> residuals;
      EPNPEstimator::Residuals(points2D, points3D, report.model, &residuals);
      for (size_t i = 0; i < residuals.size(); ++i) {
        BOOST_CHECK(residuals[i] < 1e-3);
      }

      // Test residuals of faulty points.
      EPNPEstimator::Residuals(points2D, points3D_faulty, report.model,
                               &residuals);
      for (size_t i = 0; i < residuals.size(); ++i) {
        BOOST_CHECK(residuals[i] > 0.1);
      }
    }
  }
}

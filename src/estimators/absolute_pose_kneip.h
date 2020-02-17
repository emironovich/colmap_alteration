//
// Created by elizaveta on 13.12.2019.
//

#ifndef COLMAP_ABSOLUTE_POSE_KNEIP_H
#define COLMAP_ABSOLUTE_POSE_KNEIP_H

#include "estimators/absolute_pose.h"

namespace colmap {
class P3PEstimatorKneip : public P3PEstimator {
  // Estimate the most probable solution of the P3P problem from a set of
  // three 2D-3D point correspondences.
  //
  // @param points2D   Normalized 2D image points as 3x2 matrix.
  // @param points3D   3D world points as 3x3 matrix.
  //
  // @return           Solutions for poses as a vector of a 3x4 matrix.
  //                   Number of solutions is not grater than 4.
 public:
  static std::vector<M_t> Estimate(const std::vector<X_t>& points2D,
                                   const std::vector<Y_t>& points3D);
};
}  // namespace colmap

#endif  // COLMAP_ABSOLUTE_POSE_KNEIP_H

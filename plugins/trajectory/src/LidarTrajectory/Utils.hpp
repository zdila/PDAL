#if !defined(LIDARTRAJECTORY_UTILS_HPP)
#define LIDARTRAJECTORY_UTILS_HPP 1

#include <Eigen/Dense>

namespace LidarTrajectory {

  class Utils {
  public:
    static Eigen::Matrix3d PerpProjector(const Eigen::Vector3d& v, double d) {
      // Assume v is a unit vector
      double v2 = v(0) * v(0) + v(1) * v(1);
      Eigen::Matrix3d M;
      // See proj.mac for derivation
      M(0,0) = (v(0)*v(0) * v(2) + v(1)*v(1)) / v2;
      M(1,1) = (v(1)*v(1) * v(2) + v(0)*v(0)) / v2;
      M(0,1) = -(1 - v(2)) * v(0) * v(1) / v2; M(1,0) = M(0,1);
      M(0,2) = -v(0); M(1,2) = -v(1);
      M.row(2) = v.transpose();
      // Scale first two rows to get projection error at v
      M.row(0) *= d;
      M.row(1) *= d;
      return M;
    }

  };

}
#endif

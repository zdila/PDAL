#include <GeographicLib/Math.hpp>
#include <LidarTrajectory/Pulse.hpp>

namespace LidarTrajectory {

  Pulse::Pulse(double _t, const Eigen::Vector3d& _r, double _ang)
    : t(_t)
    , r(_r)
    , n(Eigen::Vector3d::Zero())
    , d(0)
    , ang(_ang * GeographicLib::Math::degree())
  {}

  Pulse::Pulse(double _t,
               const Eigen::Vector3d& _f,
               const Eigen::Vector3d& _l,
               double _ang)
    : t(_t)
    , r((_f + _l) /2)
    , n(_f - _l)
    , d(n.norm() / 2)
    , ang(_ang * GeographicLib::Math::degree())
  {
    if (d != 0)
      n /= 2*d;
  }

}

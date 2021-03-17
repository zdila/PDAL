#if !defined(LASERTRAJ_HPP)
#define LASERTRAJ_HPP 1

#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <utility>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Core>

class position {
public:
  double t;
  Eigen::Vector3d r;
  position()
    : t(0), r(Eigen::Vector3d::Zero())
  {}
  position(double t_ = 0, double x_ = 0, double y_ = 0, double z_ = 0)
    : t(t_), r(x_, y_, z_)
  {}
  position(double t_, const Eigen::Vector3d& r_)
    : t(t_), r(r_)
  {}
};

class pulse {
public:
  double t, scan;
  Eigen::Vector3d f, l;
  pulse()
    : t(0), scan(0), f(Eigen::Vector3d::Zero()), l(Eigen::Vector3d::Zero())
  {}
  pulse(double t_,
        double xf_, double yf_, double zf_,
        double xl_, double yl_, double zl_,
        double scan_ = 0)
    : t(t_)
    , scan(scan_)
    , f(xf_, yf_, zf_)
    , l(xl_, yl_, zl_)
  {}
  pulse(double t_,
        const Eigen::Vector3d& f_, const Eigen::Vector3d& l_,
        double scan_ = 0)
    : t(t_), scan(scan_), f(f_), l(l_)
  {}
};

class trajfit {
public:
  int num;                      // Number of blocks
  double tblock, tstart;        // tend = tstart + num*tblock
  std::vector<Eigen::Vector3d> r; // positions at block boundaries
  std::vector<Eigen::Vector3d> v; // velocity*tblock at block boundaries
  // return position and velocity at time t
  std::pair<Eigen::Vector3d, Eigen::Vector3d> coords(double t) const;
  trajfit(int _num = -1, double _tblock = 1, double _tstart = 0)
    : num(_num)
    , tblock(_tblock)
    , tstart(_tstart)
    , r(num+1)                // num+1 to get both endpoints
    , v(num+1)
  {}
};

class lasertraj {
public:
  static const int ord = 4;     // For cubic
  double torg;
  Eigen::Vector3d rorg;
private:
  // fit order in N-1
  static const int N = 8;
  std::vector<double> times;
  std::vector<pulse> pulses;
  std::vector<position> esttraj;
  double dt, dtblank, dttraj, tblock, dr, dang;
  int vlevel, niter;
  Eigen::Matrix<double, N, 3> C;
  trajfit traj;                 // the fit for x,y,z
  trajfit attitude;             // the fit for yaw,pitch,roll
  template <typename T>
  static  T cubic(const T* const c, T t) {
    return  t * (t * (t * c[3] + c[2]) + c[1]) + c[0];
  }
  // Return [min, max) indexes for fit at time t.  N.B. t adjusted to be within
  // time limits of array
  template <typename T>
  static std::pair<int, int> indrange(double& t, double tblock,
                                      const std::vector<T>& array);
public:
  lasertraj(std::istream& in);
  void estimate();
  void initialfit();
  void print(std::ostream& out) const;
  Eigen::Matrix<double, ord, 3> lsfit(double t) const;
  Eigen::Matrix<double, ord, 3> ceresfit(double t, bool alt) const;
  void finalfit();
  void attitudefit();
  void attitudefitfull();
  // Projector onto plane perpendicular to v
  static Eigen::Matrix3d PerpProjector(const Eigen::Vector3d& v);
  std::pair<Eigen::Vector3d, Eigen::Vector3d> coords(double t) const {
    return traj.coords(t);
  }
  std::pair<Eigen::Vector3d, Eigen::Vector3d> att(double t) const {
    auto p = attitude.coords(t);
    double deg = std::atan2(0.0, -1.0)/180;
    p.first /= deg;
    p.second /= deg;
    return p;
  }
  template <typename T>
  static T EndPointCubic(const T& rm, const T& vm,
                         const T& rp, const T& vp,
                         const T& t, T* v = nullptr);

  static Eigen::Vector3d cubic3(const Eigen::Matrix<double, ord, 3>& C,
                                double t) {
    Eigen::Vector3d r(Eigen::Vector3d::Zero());
    for (int k = ord; k > 0;) {
      --k;
      r = t * r + C.row(k).transpose();
    }
    return r;
  }
  static Eigen::Vector3d cubic3d(const Eigen::Matrix<double, ord, 3>& C,
                                 double t) {
    Eigen::Vector3d r(Eigen::Vector3d::Zero());
    for (int k = ord; k > 1;) {
      --k;
      r = t * r + k * C.row(k).transpose();
    }
    return r;
  }
};

#endif

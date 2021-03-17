// Input created by

// las2txt --parse ptxyzanr -i xx.las -o xx.txt

// Current assumptions:

// Only first point source ID used

// Data is assumed to be ordered by GPS time

// Consider only n > 1 and r == 1 or r == n

// Sample data every dt = 0.001 s taking max height difference in each dt
// interval.

// Fit cubic path to trajectory in blocks of Dt = 0.5 s
//  x = a + b*t + c*t^2 + d*t^3

// Do the cubic path fit with three adjacent blocks

// Quantization error in positions is 0.01 m
#include "lasertraj.hpp"

#include <chrono>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <glog/logging.h>

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::NumericDiffCostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

lasertraj::lasertraj(std::istream& in)
  : torg(0)
  , rorg(Eigen::Vector3d::Zero())
  , dt(0.001)
  , dtblank(dt/20)
  , dttraj(50*dt)
  , tblock(1)
  , dr(0.01)                          // 1 cm
  , dang(1 * std::atan2(0.0, -1.0) / 180) // 1 deg
  , vlevel(0)
  , niter(50)
{
  std::string line;
  bool first = true;
  std::istringstream str;
  double told = -1;
  Eigen::Vector3d r, rf(Eigen::Vector3d::Zero());
  while (std::getline(in, line)) {
    str.clear(); str.str(line);
    int p, ret, num;
    double t, a;
    if (!(str >> p >> t
          >> r(0) >> r(1) >> r(2)
          >> a >> num >> ret))
      throw std::runtime_error("Short read on " + line);
    if (first) {
      torg = std::floor(t/10)*10;
      rorg = (Eigen::floor(r.array()/100)*100).matrix();
      first = false;
    }
    t -= torg; r -= rorg;
    times.push_back(t);
    // Only consider first and last returns
    if (! (num > 1 && (ret == 1 || ret == num)) ) continue;
    if (t < told)
      throw std::runtime_error("Time must not decrease");
    if (t == told) {
      pulse P(t, rf, r, a);
      if (pulses.size() == 0)
        pulses.push_back(P);
      else if (std::floor(t/dt) > std::floor(pulses.back().t/dt)) {
        if (t >= pulses.back().t + dtblank)
          pulses.push_back(P);
      } else {
        if (std::abs(rf(2) - r(2)) >
            std::abs(pulses.back().f(2) - pulses.back().l(2)))
          pulses.back() = P;
      }
      told = std::nextafter(t, t + 1);
    } else {
      told = t; rf = r;
    }
  }
}

void lasertraj::estimate() {
  // equation for ray
  //   (x-xl) = px * (z-zl) or x - px * z = xl - px * zl
  //   (y-yl) = py * (z-zl) or y - py * z = yl - py * zl
  //   px = (xf-xl)/(zf-zl); py = (yf-yl)/(zf-zl)
  // As a least squares problem
  //   A . [x,y,z]' = B
  // where
  //   A = [1, 0, -px]
  //       [0, 1, -py]
  //       ... for all rays
  //   B = [xl - px * zl]
  //       [yl - py * zl]
  //       ... for all rays
  int npulses = int(pulses.size());
  if (npulses == 0)
    throw std::runtime_error("No multiple returns");
  esttraj.clear();
  int j = 0;
  while (j < npulses) {
    double tstart = std::floor(pulses[j].t/dttraj) * dttraj,
      tend = tstart + dttraj;
    int jstart = j;
    while (true) {
      ++j;
      if (j == npulses || pulses[j].t >= tend)
        break;
    }
    int jend = j;
    if (jend - jstart < 2)
      continue;
    Eigen::Matrix<double, Eigen::Dynamic, 3> A(2*(jend - jstart), 3);
    Eigen::VectorXd B(2*(jend - jstart));
    double t = 0;
    for (int j = jstart, i = 0; j < jend; ++j, ++i) {
      const pulse& p = pulses[j];
      t += p.t;
      double
        px = (p.f(0) - p.l(0)) / (p.f(2) - p.l(2)),
        py = (p.f(1) - p.l(1)) / (p.f(2) - p.l(2));
      A(2*i+0, 0) = 1; A(2*i+0, 1) = 0; A(2*i+0, 2) = -px;
      A(2*i+1, 0) = 0; A(2*i+1, 1) = 1; A(2*i+1, 2) = -py;
      B(2*i+0) = p.l(0) - px * p.l(2);
      B(2*i+1) = p.l(1) - py * p.l(2);
    }
    t /= (jend - jstart);
    Eigen::Vector3d pos(A.jacobiSvd(Eigen::ComputeThinU |
                                    Eigen::ComputeThinV).solve(B));
    esttraj.push_back(position(t, pos));
  }
}

void lasertraj::initialfit() {
  // Fit cubic to esttraj
  // x = ax + bx * t + cx * t^2 + dx * t^3, etc
  // 3 indep least squares problem (for x y z)
  //   A . [ax, bx, cx, dx] = B
  // where
  //   A = [1, t, t^2, t^3]
  //       ... for all times
  //   B = [x]
  //       ... for all times

  int num = int(esttraj.size());
  Eigen::Matrix<double, Eigen::Dynamic, N> A(num, N);
  Eigen::Matrix<double, Eigen::Dynamic, 3> B(num, 3);
  int j = 0;
  for (auto p = esttraj.begin(); p < esttraj.end(); ++p, ++j) {
    double t = p->t;
    A(j, 0) = 1;
    for (int k = 1; k < N; ++k)
      A(j, k) = t * A(j, k-1);
    B.row(j) = p->r.transpose();
  }
  C = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B);
}

void lasertraj::print(std::ostream& out) const {
  out << std::fixed;
  for (auto p = esttraj.begin(); p < esttraj.end(); ++p) {
    double t = p->t;
    Eigen::Vector3d s = p->r + rorg;
    Eigen::Vector3d r(Eigen::Vector3d::Zero());
    for (int k = N; k > 0;) {
      --k;
      r = t * r + C.row(k).transpose();
    }
    r = r + rorg;
    out << std::setprecision(4) << t + torg << " "
        << std::setprecision(3)
        << s(0) << " " << s(1) << " " << s(2) << " "
        << r(0) << " " << r(1) << " " << r(2) << "\n";
  }
}

  // Return min/max indexes for fit at time t.  N.B. t adjusted to be within
  // time limits of esttraj
template <typename T>
std::pair<int, int> lasertraj::indrange(double& t, double tblock,
                                        const std::vector<T>& array) {
  int ntraj = int(array.size());
  if (ntraj == 0)
    throw std::runtime_error("No trajectory points");
  t = std::max(array[0].t, t);
  t = std::min(array[ntraj-1].t, t);
  double tmin = t - tblock, tmax = t + tblock;
  int jstart = ntraj, jend = -1;
  int j = 0;
  for (auto p = array.begin(); p < array.end(); ++p, ++j) {
    if (p->t > tmin && p->t < tmax) {
      jstart = std::min(j, jstart);
      jend = std::max(j, jend);
    }
  }
  ++jend;
  return std::make_pair(jstart, jend);
}

Eigen::Matrix<double, lasertraj::ord, 3>
lasertraj::lsfit(double t) const {
  std::pair<int, int> range = indrange(t, tblock, esttraj);
  int jstart = range.first, jend = range.second, num = jend - jstart;
  Eigen::Matrix<double, Eigen::Dynamic, ord> A(num, ord);
  Eigen::Matrix<double, Eigen::Dynamic, 3> B(num, 3);
  for (int i = 0, j = jstart; j < jend; ++i, ++j) {
    const position& p = esttraj[j];
    double t = p.t;
    A(i, 0) = 1;
    for (int k = 1; k < ord; ++k)
      A(i, k) = t * A(i, k-1);
    B.row(i) = p.r.transpose();
  }
  return A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B);
}

class RayError {
  const pulse _p;
public:
  RayError(const pulse& p) : _p(p) {}
  template <typename T>
  bool operator()(const T* const traj, // ord * 3 comp cubic fit to x,y,z
                  const T* const slope, // inverse slopes for x,y
                  const T* const range, // zranges for first,last
                  // 6 residuals
                  T* residuals) const {
    const T t = T(_p.t);
    // Platform position
    T x = t * (t * (t * traj[0+3] + traj[0+2]) + traj[0+1]) + traj[0+0];
    T y = t * (t * (t * traj[4+3] + traj[4+2]) + traj[4+1]) + traj[4+0];
    T z = t * (t * (t * traj[8+3] + traj[8+2]) + traj[8+1]) + traj[8+0];
    const T f[3] = { T(_p.f(0)), T(_p.f(1)), T(_p.f(2)) };
    const T l[3] = { T(_p.l(0)), T(_p.l(1)), T(_p.l(2)) };
    T
      fx = x + slope[0] * range[0],
      fy = y + slope[1] * range[0],
      fz = z +            range[0],
      lx = x + slope[0] * range[1],
      ly = y + slope[1] * range[1],
      lz = z +            range[1];
    residuals[0] = fx - f[0];
    residuals[1] = fy - f[1];
    residuals[2] = fz - f[2];
    residuals[3] = lx - l[0];
    residuals[4] = ly - l[1];
    residuals[5] = lz - l[2];
    return true;
  }
};

// Alternate version of RayError which does not require slopes and ranges to be
// estimated.  Instead use the fact that the "best" ray joins the platform to
// the mean of the first and last returns.
class RayErrorAlt {
  double _t;
  const Eigen::Vector3d _mean, _diff;
public:
  RayErrorAlt(const pulse& p)
    : _t(p.t)
    , _mean((p.l + p.f)/2)
    , _diff((p.l - p.f)/2)
 {}
  template <typename T>
  bool operator()(const T* const traj, // ord * 3 comp cubic fit to x,y,z
                  // 3 residuals
                  T* residuals) const {
    using std::sqrt;
    const T t = T(_t),
      xm = T(_mean(0)), ym = T(_mean(1)), zm = T(_mean(2)),
      xd = T(_diff(0)), yd = T(_diff(1)), zd = T(_diff(2));
    // Platform position
    T x = t * (t * (t * traj[0+3] + traj[0+2]) + traj[0+1]) + traj[0+0];
    T y = t * (t * (t * traj[4+3] + traj[4+2]) + traj[4+1]) + traj[4+0];
    T z = t * (t * (t * traj[8+3] + traj[8+2]) + traj[8+1]) + traj[8+0];
    x -= xm; y -= ym; z -= zm;  // relative to mean return
    T d = sqrt( x*x + y*y + z*z );
    x /= d; y /= d; z /= d;     // convert to unit vector
    // cross product of [x,y,z] and [xd,yd,zd]
    residuals[0] = y * zd - z * yd;
    residuals[1] = z * xd - x * zd;
    residuals[2] = x * yd - y * xd;
    return true;
  }
};

template <typename T>
T lasertraj::EndPointCubic(const T& rm, const T& vm,
                           const T& rp, const T& vp,
                           const T& t, T* v) {
  T rs = rp + rm,
    rd = rp - rm,
    vs = vp + vm,
    vd = vp - vm,
    a0 = (4.0 * rs - vd) / 8.0,
    a1 = (6.0 * rd - vs) / 4.0,
    a2 =        vd       / 2.0,
    a3 = -2.0 * rd + vs       ;
  if (v)
    *v = t * (t * 3.0 * a3 + 2.0 * a2) + a1;
  return t * (t * (t * a3 + a2) + a1) + a0;
}

Eigen::Matrix3d lasertraj::PerpProjector(const Eigen::Vector3d& v) {
  double va = v.norm();
  Eigen::Vector3d vn = v/va;
  double v2 = vn(0) * vn(0) + vn(1) * vn(1);
  Eigen::Matrix3d M;
  // See proj.mac for derivation
  M(0,0) = (vn(0)*vn(0) * vn(2) + vn(1)*vn(1)) / v2;
  M(1,1) = (vn(1)*vn(1) * vn(2) + vn(0)*vn(0)) / v2;
  M(0,1) = -(1 - vn(2)) * vn(0) * vn(1) / v2; M(1,0) = M(0,1);
  M(0,2) = -vn(0); M(1,2) = -vn(1);
  M.row(2) = vn.transpose();
    // Scale first two rows to get projection error at v
  M.row(0) *= va;
  M.row(1) *= va;
  return M;
}
// Modify (and simplify) by projecting platform
class RayErrorAlt2 {
  double _t;
  const Eigen::Vector3d _mean;
  Eigen::Matrix3d _M;
public:
  RayErrorAlt2(const pulse& p)
    : _t(p.t)
    , _mean((p.l + p.f)/2)
  {
    double
      px = p.f(0) - p.l(0), py = p.f(1) - p.l(1), pz = p.f(2) - p.l(2),
      pa = std::sqrt(px*px + py*py + pz*pz);
    px /= pa; pz /= pa; py /= pa;
    double p2 = px*px + py*py;
    // See proj.mac for derivation
    _M(0,0) = (px*px * pz + py*py) / p2; _M(0,1) = -(1 - pz) * px * py / p2;
    _M(1,1) = (px*px + py*py * pz) / p2; _M(1,0) = _M(0,1);
    _M(0,2) = -px; _M(1,2) = -py; _M(2,0) = px; _M(2,1) = py; _M(2,2) = pz;
    // Scale first two rows to get projection error at first return
    pa /= 2;
    _M(0,0) *= pa; _M(0,1) *= pa; _M(0,2) *= pa;
    _M(1,0) *= pa; _M(1,1) *= pa; _M(1,2) *= pa;
  }
  template <typename T>
  bool operator()(const T* const traj, // ord * 3 comp cubic fit to x,y,z
                  // 2 residuals
                  T* residuals) const {
    const T t = T(_t);
    // Platform position
    T x = t * (t * (t * traj[0+3] + traj[0+2]) + traj[0+1]) + traj[0+0];
    T y = t * (t * (t * traj[4+3] + traj[4+2]) + traj[4+1]) + traj[4+0];
    T z = t * (t * (t * traj[8+3] + traj[8+2]) + traj[8+1]) + traj[8+0];
    // relative to mean return
    x -= T(_mean(0)); y -= T(_mean(1)); z -= T(_mean(2));
    // Transform
    T xt = T(_M(0,0)) * x + T(_M(0,1)) * y + T(_M(0,2)) * z;
    T yt = T(_M(1,0)) * x + T(_M(1,1)) * y + T(_M(1,2)) * z;
    T zt = T(_M(2,0)) * x + T(_M(2,1)) * y + T(_M(2,2)) * z;
    // Project
    residuals[0] = xt / zt;
    residuals[1] = yt / zt;
    return true;
  }
};

Eigen::Matrix<double, lasertraj::ord, 3> lasertraj::ceresfit(double t,
                                                             bool alt) const {
  ceres::Problem problem;
  std::pair<int, int> range = indrange(t, tblock, pulses);
    int jstart = range.first, jend = range.second, num = jend - jstart;
  // Unknowns are
  // cubic fit for x,y,z
  Eigen::Matrix<double, ord, 3> C = lsfit(t);
  std::chrono::steady_clock::time_point start =
    std::chrono::steady_clock::now();
  if (vlevel > 0) {
    std::cout << "Number of pulses " << num << "\n";
    std::cout << "Initial fit\n" << C << "\n";
    std::cout << "Initial pos + velocity\n"
              << std::fixed << std::setprecision(3)
              << (rorg + cubic3(C, t)).transpose() << "\n"
              << cubic3d(C, t).transpose() << "\n";
  }
  int skip = 1;
  num /= skip;
  std::vector<double> slope(alt ? 0 : 2*num), // x,y inverse slopes
    zrange(alt ? 0 : 2*num);                  // first,last ranges in z
  if (!alt) {
    // Fill in initial estimates for slope and range
    for (int i = 0; i < num; ++i) {
      const pulse& p = pulses[i*skip + jstart];
      double t0 = p.t;
      // est platform position
      Eigen::Vector3d r0 = cubic3(C, t0);
      // ray to midpoint of returns
      Eigen::Vector3d rm = (p.f + p.l) / 2 - r0;
      slope[2 * i + 0] = rm(0) / rm(2);
      slope[2 * i + 1] = rm(1) / rm(2);
      zrange[2 * i + 0] = p.f(2) - r0(2);
      zrange[2 * i + 1] = p.l(2) - r0(2);
    }
  }
  const double pi = std::atan2(0.0, -1.0);
  // Set up residual blocks
  for (int i = 0; i < num; ++i) {
    const pulse& p = pulses[i*skip + jstart];
    double t0 = p.t;
    // Apply full cosine window to [t-tblock, t+tblock]
    double weight = (1 + std::cos((t0 - t) / tblock * pi)) / 2;
    if (!alt) {
      ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<RayError,
                                        6,     // number of residuals
                                        ord*3, // components in C
                                        2, 2>  // slopes + ranges
        (new RayError(p));
      ceres::LossFunction* loss_function =
        (ceres::LossFunction*)
        (new ceres::ScaledLoss(new ceres::CauchyLoss(dr), weight,
                               ceres::TAKE_OWNERSHIP));

      problem.AddResidualBlock(cost_function, loss_function,
                               C.data(), // laid out column major order
                               slope.data() + 2*i,
                               zrange.data() + 2*i);
    } else {
      ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<RayErrorAlt,
                                        3,     // number of residuals
                                        ord*3> // components in C
        (new RayErrorAlt(p));
      ceres::LossFunction* loss_function =
        (ceres::LossFunction*)
        (new ceres::ScaledLoss(new ceres::CauchyLoss(dr), weight,
                               ceres::TAKE_OWNERSHIP));

      problem.AddResidualBlock(cost_function, loss_function,
                               C.data()); // laid out column major order
    }
  }
  if (vlevel > 0) {
    int nresid = (alt ? 3 : 6) * num;
    // total number of residuals
    std::vector<double> resid(nresid);
    ceres::Problem::EvaluateOptions opt;
    opt.apply_loss_function = false;
    double cost;
    problem.Evaluate(opt, &cost, &resid, NULL, NULL);
    double sum = 0;
    for (int i = 0; i < nresid; ++i) sum += resid[i] * resid[i];
    std::cerr << "Initial cost " << cost
              << " RMS resid " << std::sqrt(sum/nresid) << "\n";
  }

  // Make Ceres automatically detect the bundle structure. Note that the
  // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is
  // slower for standard bundle adjustment problems.
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_SCHUR;
  if (vlevel > 0)
    options.minimizer_progress_to_stdout = true;
  else
    options.logging_type = ceres::SILENT;
  options.max_linear_solver_iterations = niter;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  if (vlevel > 0)
    std::cout << summary.FullReport() << "\n";

  if (vlevel > 0) {
    int nresid = (alt ? 3 : 6) * num;
    // total number of residuals
    std::vector<double> resid(nresid);
    ceres::Problem::EvaluateOptions opt;
    opt.apply_loss_function = false;
    double cost;
    problem.Evaluate(opt, &cost, &resid, NULL, NULL);
    double sum = 0;
    for (int i = 0; i < nresid; ++i) sum += resid[i] * resid[i];
    std::cerr << "Initial cost " << cost
              << " RMS resid " << std::sqrt(sum/nresid) << "\n";
  }

  if (vlevel > 0) {
  std::cout << "Final fit\n" << C << "\n";
  std::cout << "Final pos + velocity\n"
            << std::fixed << std::setprecision(3)
            << (rorg + cubic3(C, t)).transpose() << "\n"
            << cubic3d(C, t).transpose() << "\n";
  }
  std::chrono::steady_clock::time_point end =
    std::chrono::steady_clock::now();
  std::cerr << "elapsed time "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms\n";
  return C;
}

//  RayErrorAlt adapted so that cubic is speficfied by values and derivatives
// at end points;
class RayErrorEndpoint {
  double _t;
  const Eigen::Vector3d _mean, _diff;
public:
  RayErrorEndpoint(const pulse& p, double t)
    : _t(t)                     // fractional time in block
    , _mean((p.l + p.f)/2)
    , _diff((p.l - p.f)/2)
 {}
  template <typename T>
  bool operator()(const T* const rm, // 3 vec for pos at beg
                  const T* const vm, // 3 vec for vel at beg
                  const T* const rp, // 3 vec for pos at end
                  const T* const vp, // 3 vec for vel at end
                  // 3 residuals
                  T* residuals) const {
    using std::sqrt;
    const T t = T(_t),
      m[3] = { T(_mean(0)), T(_mean(1)), T(_mean(2)) },
      d[3] = { T(_diff(0)), T(_diff(1)), T(_diff(2)) };
    T p[3];                     // Platform direction
    for (int i = 0; i < 3; ++i) {
      T rs = rp[i] + rm[i],
        rd = rp[i] - rm[i],
        vs = vp[i] + vm[i],
        vd = vp[i] - vm[i],
        a0 = (4.0 * rs - vd) / 8.0,
        a1 = (6.0 * rd - vs) / 4.0,
        a2 =        vd       / 2.0,
        a3 = -2.0 * rd + vs       ;
      p[i] = t * (t * (t * a3 + a2) + a1) + a0 - m[i];

    }
    T s = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
    p[0] /= s; p[1] /= s; p[2] /= s;     // convert to unit vector
    // The "right" formula for the residue is R1 = (I - pp).d; where (I - pp)
    // is the perpendicular projector.  Instead we compute the cross product R2
    // = p x d.  Since R1 = R2 x p, this makes no difference as far as the
    // residue minimization goes.
    residuals[0] = p[1] * d[2] - p[2] * d[1];
    residuals[1] = p[2] * d[0] - p[0] * d[2];
    residuals[2] = p[0] * d[1] - p[1] * d[0];
    return true;
  }
};

//  RayErrorEndpoint using projection to 2d plane
class RayErrorEndpoint2 {
  double _t;
  const Eigen::Vector3d _mean;
  Eigen::Matrix3d _M;
public:
  RayErrorEndpoint2(const pulse& p, double t)
    : _t(t)                     // fractional time in block
    , _mean((p.f + p.l)/2)
    , _M(lasertraj::PerpProjector((p.f - p.l)/2))
  {}
  template <typename T>
  bool operator()(const T* const rm, // 3 vec for pos at beg
                  const T* const vm, // 3 vec for vel at beg
                  const T* const rp, // 3 vec for pos at end
                  const T* const vp, // 3 vec for vel at end
                  // 2 residuals
                  T* residuals) const {
    const T t = T(_t),
      m[3] = { T(_mean(0)), T(_mean(1)), T(_mean(2)) };
    T p[3];                     // Platform direction
    for (int i = 0; i < 3; ++i)
      p[i] = lasertraj::EndPointCubic(rm[i], vm[i], rp[i], vp[i], t) - m[i];
    // Transform
    T xt = T(_M(0,0)) * p[0] + T(_M(0,1)) * p[1] + T(_M(0,2)) * p[2];
    T yt = T(_M(1,0)) * p[0] + T(_M(1,1)) * p[1] + T(_M(1,2)) * p[2];
    T zt = T(_M(2,0)) * p[0] + T(_M(2,1)) * p[1] + T(_M(2,2)) * p[2];
    // Project
    residuals[0] = xt / zt;
    residuals[1] = yt / zt;
    return true;
  }
};

class AccelJumpConstraint {
public:
  AccelJumpConstraint() {}
  template <typename T>
  bool operator()(const T* const ra, // 3 vec for pos at beg
                  const T* const va, // 3 vec for vel at beg
                  const T* const vb, // 3 vec for vel at cent
                  const T* const rc, // 3 vec for pos at end
                  const T* const vc, // 3 vec for vel at end
                  // 3 residuals
                  T* residual) const {
    // For no jump in the acceleration between a-b and b-c fits we want
    // vb = (3*(rc-ra) - (vc+va)) / 4
    for (int i = 0; i < 3; ++i)
      residual[i] = vb[i] - (3.0 * (rc[i] - ra[i]) - (vc[i] + va[i])) / 4.0;
    return true;
  }
};

void lasertraj::finalfit() {
  ceres::Problem problem;
  std::chrono::steady_clock::time_point start =
    std::chrono::steady_clock::now();

  int npulses = int(pulses.size());
  if (npulses == 0)
    throw std::runtime_error("No multiple returns");
  traj.tblock = tblock;
  double tstart = std::floor(pulses[0].t/tblock) * tblock;
  traj.tstart = tstart;
  int num = std::ceil(pulses[npulses - 1].t/tblock)
    - std::floor(pulses[0].t/tblock) - 1;
  traj.num = num;
  traj.r.resize(num+1);
  traj.v.resize(num+1);
  for (int i = 0; i <= num; ++i) {
    double t = tstart + i * tblock;
    Eigen::Matrix<double, ord, 3> C = lsfit(t);
    traj.r[i] = cubic3(C, t);
    traj.v[i] = cubic3d(C, t) * tblock;
  }
  int skip = 1;
  int k = npulses / skip;
  // Set up residual blocks for pulses
  for (int j = 0; j < k; ++j) {
    const pulse& p = pulses[j*skip];
    double t = p.t;
    int i = std::min(num-1,
                     std::max(0, int(std::floor((t - tstart) / tblock))));
    double tf = (t - tstart) / tblock - (i + 0.5);
    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<RayErrorEndpoint2,
                                      2,       // number of residuals
                                      3,3,3,3> // data for cubic fit
      (new RayErrorEndpoint2(p, tf));
    ceres::LossFunction* loss_function =
      (ceres::LossFunction*)(new ceres::CauchyLoss(dr));
    problem.AddResidualBlock(cost_function, loss_function,
                             traj.r[i  ].data(),
                             traj.v[i  ].data(),
                             traj.r[i+1].data(),
                             traj.v[i+1].data());
  }
  // Not sure what a reasonable value is here
  double acceljumpweight = 1;
  for (int i = 1; i < num; ++i) {
    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<AccelJumpConstraint,
                                      3,       // number of residuals
                                      3,3,3,3,3> // data for cubic fit
      (new AccelJumpConstraint());
    ceres::LossFunction* loss_function =
      (ceres::LossFunction*)
      (new ceres::ScaledLoss(NULL, acceljumpweight,
                             ceres::TAKE_OWNERSHIP));
    problem.AddResidualBlock(cost_function, loss_function,
                             traj.r[i-1].data(),
                             traj.v[i-1].data(),
                             traj.v[i  ].data(),
                             traj.r[i+1].data(),
                             traj.v[i+1].data());
  }
  if (vlevel > 0) {
   int nresid = 3 * (k + num-1);
    // total number of residuals
    std::vector<double> resid(nresid);
    ceres::Problem::EvaluateOptions opt;
    opt.apply_loss_function = false;
    double cost;
    problem.Evaluate(opt, &cost, &resid, NULL, NULL);
    double sum = 0;
    for (int i = 0; i < nresid; ++i) sum += resid[i] * resid[i];
    std::cerr << "Initial cost " << cost
              << " RMS resid " << std::sqrt(sum/nresid) << "\n";
  }

  // Make Ceres automatically detect the bundle structure. Note that the
  // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is
  // slower for standard bundle adjustment problems.
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_SCHUR;
  if (vlevel > 0)
    options.minimizer_progress_to_stdout = true;
  else
    options.logging_type = ceres::SILENT;
  options.max_linear_solver_iterations = niter;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  if (vlevel > 0)
    std::cout << summary.FullReport() << "\n";

  if (vlevel > 0) {
    int nresid = 3 * (k + num - 1);
    // total number of residuals
    std::vector<double> resid(nresid);
    ceres::Problem::EvaluateOptions opt;
    opt.apply_loss_function = false;
    double cost;
    problem.Evaluate(opt, &cost, &resid, NULL, NULL);
    double sum = 0;
    for (int i = 0; i < nresid; ++i) sum += resid[i] * resid[i];
    std::cerr << "Initial cost " << cost
              << " RMS resid " << std::sqrt(sum/nresid) << "\n";
  }

  std::chrono::steady_clock::time_point end =
    std::chrono::steady_clock::now();
  std::cerr << "elapsed time "
            << std::chrono::duration_cast<std::chrono::milliseconds>
    (end - start).count()
            << "ms\n";

}

//  Attitude error calculated assuming the position is given
class AttitudeError {
  double _scan, _t;
  const Eigen::Vector3d _dir;
public:
  AttitudeError(const pulse& p,
                double t, const Eigen::Vector3d r)
    : _scan(p.scan * std::atan2(0.0, -1.0) / 180)
    , _t(t)                     // fractional time in block
    , _dir(((p.f + p.l)/2 - r).normalized())
  {}
  template <typename T>
  bool operator()(const T* const am, // 3 vec for attitude  at beg
                  const T* const bm, // 3 vec for attitude' at beg
                  const T* const ap, // 3 vec for attitude  at end
                  const T* const bp, // 3 vec for attitude' at end
                  // 2 residuals
                  T* residuals) const {
    using std::sin; using std::cos;
    // Only use first two components (yaw, pitch) of attitude roll is
    // effectively zero (definition of scan angle includes the roll of the
    // platform).
    const int ncomps = 2;
    const T t = T(_t), scan = T(_scan);
    T a[ncomps];                     // Platform attitude
    for (int i = 0; i < ncomps; ++i)
      a[i] = lasertraj::EndPointCubic(am[i], bm[i], ap[i], bp[i], t);
    T d0[3] = { T(_dir(0)), T(_dir(1)), T(_dir(2)) };
    // Apply rotz(yaw)
    T sz = sin(a[0]), cz = cos(a[0]);
    T d1[3] = { cz * d0[0] - sz * d0[1],
                sz * d0[0] + cz * d0[1],
                d0[2] };
    // Apply rotx(-pitch)
    T sx = -sin(a[1]), cx = cos(a[1]);
    T d2[3] = { d1[0],
                cx * d1[1] - sx * d1[2],
                sx * d1[1] + cx * d1[2] };
    // Apply roty(scan-roll)
    // scan sign is wrong; it seems to be positive on left.  scan includes the
    // roll of the platform so there's no need to include a[2].  So...
    // apply roty(-scan)
    T sy = -sin(scan), cy = cos(scan);
    residuals[0] = sy * d2[2] + cy * d2[0];
    residuals[1] = d2[1];
    return true;
  }
};

void lasertraj::attitudefit() {
  // Assume finalfit has initialize traj
  ceres::Problem problem;
  std::chrono::steady_clock::time_point start =
    std::chrono::steady_clock::now();

  int npulses = int(pulses.size());
  if (npulses == 0)
    throw std::runtime_error("No multiple returns");
  attitude = traj;
  double deg = std::atan2(0.0, -1.0) / 180;
  for (int j = 0; j <= traj.num; ++j) {
    attitude.r[j] = Eigen::Vector3d(-90 * deg, 0.0, 0.0);
    attitude.v[j] = Eigen::Vector3d::Zero();
  }
  int skip = 1;
  int k = npulses / skip;
  // Set up residual blocks for pulses
  for (int j = 0; j < k; ++j) {
    const pulse& p = pulses[j*skip];
    double t = p.t;
    int i = std::min(attitude.num-1,
                     std::max(0, int(std::floor((t - attitude.tstart) /
                                                attitude.tblock))));
    double tf = (t - attitude.tstart) / attitude.tblock - (i + 0.5);
    Eigen::Vector3d r = traj.coords(t).first;
    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<AttitudeError,
                                      2,       // number of residuals
                                      3,3,3,3> // data for cubic fit
      (new AttitudeError(p, tf, r));
    ceres::LossFunction* loss_function =
      (ceres::LossFunction*)(new ceres::CauchyLoss(dang));
    problem.AddResidualBlock(cost_function, loss_function,
                             attitude.r[i  ].data(),
                             attitude.v[i  ].data(),
                             attitude.r[i+1].data(),
                             attitude.v[i+1].data());
  }
  // Not sure what a reasonable value is here
  double acceljumpweight = 1;
  for (int i = 1; i < attitude.num; ++i) {
    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<AccelJumpConstraint,
                                      3,       // number of residuals
                                      3,3,3,3,3> // data for cubic fit
      (new AccelJumpConstraint());
    ceres::LossFunction* loss_function =
      (ceres::LossFunction*)
      (new ceres::ScaledLoss(NULL, acceljumpweight,
                             ceres::TAKE_OWNERSHIP));
    problem.AddResidualBlock(cost_function, loss_function,
                             attitude.r[i-1].data(),
                             attitude.v[i-1].data(),
                             attitude.v[i  ].data(),
                             attitude.r[i+1].data(),
                             attitude.v[i+1].data());
  }
  if (vlevel > 0) {
   int nresid = 3 * (k + attitude.num-1);
    // total number of residuals
    std::vector<double> resid(nresid);
    ceres::Problem::EvaluateOptions opt;
    opt.apply_loss_function = false;
    double cost;
    problem.Evaluate(opt, &cost, &resid, NULL, NULL);
    double sum = 0;
    for (int i = 0; i < nresid; ++i) sum += resid[i] * resid[i];
    std::cerr << "Initial cost " << cost
              << " RMS resid " << std::sqrt(sum/nresid) << "\n";
  }

  // Make Ceres automatically detect the bundle structure. Note that the
  // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is
  // slower for standard bundle adjustment problems.
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_SCHUR;
  if (vlevel > 0)
    options.minimizer_progress_to_stdout = true;
  else
    options.logging_type = ceres::SILENT;
  options.max_linear_solver_iterations = niter;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  if (vlevel > 0)
    std::cout << summary.FullReport() << "\n";

  if (vlevel > 0) {
    int nresid = 3 * (k + attitude.num - 1);
    // total number of residuals
    std::vector<double> resid(nresid);
    ceres::Problem::EvaluateOptions opt;
    opt.apply_loss_function = false;
    double cost;
    problem.Evaluate(opt, &cost, &resid, NULL, NULL);
    double sum = 0;
    for (int i = 0; i < nresid; ++i) sum += resid[i] * resid[i];
    std::cerr << "Initial cost " << cost
              << " RMS resid " << std::sqrt(sum/nresid) << "\n";
  }

  std::chrono::steady_clock::time_point end =
    std::chrono::steady_clock::now();
  std::cerr << "elapsed time "
            << std::chrono::duration_cast<std::chrono::milliseconds>
    (end - start).count()
            << "ms\n";

}

//  Attitude error with position optimization
class AttitudeErrorFull {
  double _scan, _t;
  const Eigen::Vector3d _pos;
public:
  AttitudeErrorFull(const pulse& p, double t)
    : _scan(p.scan * std::atan2(0.0, -1.0) / 180)
    , _t(t)                     // fractional time in block
    , _pos((p.f + p.l)/2)
  {}
  template <typename T>
  bool operator()(const T* const rm, // 3 vec for pos at beg
                  const T* const vm, // 3 vec for vel at beg
                  const T* const rp, // 3 vec for pos at end
                  const T* const vp, // 3 vec for vel at end
                  const T* const am, // 3 vec for attitude  at beg
                  const T* const bm, // 3 vec for attitude' at beg
                  const T* const ap, // 3 vec for attitude  at end
                  const T* const bp, // 3 vec for attitude' at end
                  // 2 residuals
                  T* residuals) const {
    using std::sin; using std::cos;
    // Only use first two components (yaw, pitch) of attitude roll is
    // effectively zero (definition of scan angle includes the roll of the
    // platform).
    const int ncomps = 2;
    const T t = T(_t), scan = T(_scan);
    T d0[3],                     // Platform direction

      m[3] = { T(_pos(0)), T(_pos(1)), T(_pos(2)) };

    for (int i = 0; i < 3; ++i)
      d0[i] = lasertraj::EndPointCubic(rm[i], vm[i], rp[i], vp[i], t) - m[i];
    // Don't normalize; instead do a projection at the end
    T a[ncomps];                     // Platform attitude
    for (int i = 0; i < ncomps; ++i)
      a[i] = lasertraj::EndPointCubic(am[i], bm[i], ap[i], bp[i], t);
    // Apply rotz(yaw)
    T sz = sin(a[0]), cz = cos(a[0]);
    T d1[3] = { cz * d0[0] - sz * d0[1],
                sz * d0[0] + cz * d0[1],
                d0[2] };
    // Apply rotx(-pitch)
    //    T sx = -sin(a[1]), cx = cos(a[1]);
    T sx = T(0.0), cx = T(1.0);
    T d2[3] = { d1[0],
                cx * d1[1] - sx * d1[2],
                sx * d1[1] + cx * d1[2] };
    // Apply roty(scan-roll)
    // scan sign is wrong; it seems to be positive on left.  scan includes the
    // roll of the platform so there's no need to include a[2].  So...
    // apply roty(-scan)
    T sy = -sin(scan), cy = cos(scan);
    T d3[3] = { sy * d2[2] + cy * d2[0],
                d2[1],
                cy * d2[2] - sy * d2[0] };
    residuals[0] = d3[0] / d3[2];
    residuals[1] = d3[1] / d3[2] * 4.0;
    return true;
  }
};

//  Damp pitch
class DampPitch {
  double _targetpitch;
public:
  DampPitch(double targetpitch)
    : _targetpitch(targetpitch * std::atan2(0.0, -1.0) / 180)
  {}
  template <typename T>
  bool operator()(const T* const am, // attitude  at beg
                  const T* const bm, // attitude' at beg
                  const T* const ap, // attitude  at end
                  const T* const bp, // attitude' at end
                  // 5 residuals
                  T* residuals) const {
    for (int i = 0; i < 5; ++i) {
      T t = T(0.2 * i - 0.4);   // [-0.4, -0.2, 0.0, 0.2, 0.4];
      // only apply to pitch component (index 1)
      T p = lasertraj::EndPointCubic(am[1], bm[1], ap[1], bp[1], t);
      residuals[i] = p - _targetpitch;
    }
    return true;
  }
};

void lasertraj::attitudefitfull() {
  // Assume finalfit has initialize traj
  ceres::Problem problem;
  std::chrono::steady_clock::time_point start =
    std::chrono::steady_clock::now();

  int npulses = int(pulses.size());
  if (npulses == 0)
    throw std::runtime_error("No multiple returns");
  /*  for (int j = 0; j <= traj.num; ++j) {
      if (j > 0 && j < traj.num)
      traj.r[j] = traj.r[0] + (traj.r[traj.num] - traj.r[0]) *
      (double(j)/traj.num);
      traj.v[j] = (traj.r[traj.num] - traj.r[0]) / traj.num;
      }
  */
  attitude = traj;
  double deg = std::atan2(0.0, -1.0) / 180;
  for (int j = 0; j <= attitude.num; ++j) {
    attitude.r[j] = Eigen::Vector3d(-90 * deg, 0.0, 0.0);
    attitude.v[j] = Eigen::Vector3d::Zero();
  }
  int skip = 1;
  int k = npulses / skip;
  // Set up residual blocks for pulses
  for (int j = 0; j < k; ++j) {
    const pulse& p = pulses[j*skip];
    double t = p.t;
    int i = std::min(attitude.num-1,
                     std::max(0, int(std::floor((t - attitude.tstart) /
                                                attitude.tblock))));
    double tf = (t - attitude.tstart) / attitude.tblock - (i + 0.5);
    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<AttitudeErrorFull,
                                      2,       // number of residuals
                                      3,3,3,3, // data for cubic fit to pos
                                      3,3,3,3> // data for fit to attitude
      (new AttitudeErrorFull(p, tf));
    ceres::LossFunction* loss_function =
      (ceres::LossFunction*)(new ceres::CauchyLoss(dang));
    problem.AddResidualBlock(cost_function, loss_function,
                             traj.r[i  ].data(),
                             traj.v[i  ].data(),
                             traj.r[i+1].data(),
                             traj.v[i+1].data(),
                             attitude.r[i  ].data(),
                             attitude.v[i  ].data(),
                             attitude.r[i+1].data(),
                             attitude.v[i+1].data());
  }
  // Not sure what a reasonable value is here
  double acceljumpweight = 1;
  for (int i = 1; i < traj.num; ++i) {
    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<AccelJumpConstraint,
                                      3,       // number of residuals
                                      3,3,3,3,3> // data for cubic fit
      (new AccelJumpConstraint());
    ceres::LossFunction* loss_function =
      (ceres::LossFunction*)
      (new ceres::ScaledLoss(NULL, acceljumpweight,
                             ceres::TAKE_OWNERSHIP));
    problem.AddResidualBlock(cost_function, loss_function,
                             traj.r[i-1].data(),
                             traj.v[i-1].data(),
                             traj.v[i  ].data(),
                             traj.r[i+1].data(),
                             traj.v[i+1].data());
  }
  double attitudejumpweight = 1;
  for (int i = 1; i < attitude.num; ++i) {
    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<AccelJumpConstraint,
                                      3,       // number of residuals
                                      3,3,3,3,3> // data for cubic fit
      (new AccelJumpConstraint());
    ceres::LossFunction* loss_function =
      (ceres::LossFunction*)
      (new ceres::ScaledLoss(NULL, attitudejumpweight,
                             ceres::TAKE_OWNERSHIP));
    problem.AddResidualBlock(cost_function, loss_function,
                             attitude.r[i-1].data(),
                             attitude.v[i-1].data(),
                             attitude.v[i  ].data(),
                             attitude.r[i+1].data(),
                             attitude.v[i+1].data());
  }

  if (false) {
    double pitchdampweight = 0.05;
    double targetpitch = 2;
    for (int i = 0; i < attitude.num; ++i) {
      ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<DampPitch,
                                        5,       // number of residuals
                                        3,3,3,3> // data for cubic fit
        (new DampPitch(targetpitch));
      ceres::LossFunction* loss_function =
        (ceres::LossFunction*)
        (new ceres::ScaledLoss(NULL, pitchdampweight,
                               ceres::TAKE_OWNERSHIP));
      problem.AddResidualBlock(cost_function, loss_function,
                               attitude.r[i  ].data(),
                               attitude.v[i  ].data(),
                               attitude.r[i+1].data(),
                               attitude.v[i+1].data());
    }
  }

  if (vlevel > 0) {
   int nresid = 3 * (k + attitude.num-1);
    // total number of residuals
    std::vector<double> resid(nresid);
    ceres::Problem::EvaluateOptions opt;
    opt.apply_loss_function = false;
    double cost;
    problem.Evaluate(opt, &cost, &resid, NULL, NULL);
    double sum = 0;
    for (int i = 0; i < nresid; ++i) sum += resid[i] * resid[i];
    std::cerr << "Initial cost " << cost
              << " RMS resid " << std::sqrt(sum/nresid) << "\n";
  }

  // Make Ceres automatically detect the bundle structure. Note that the
  // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is
  // slower for standard bundle adjustment problems.
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_SCHUR;
  if (vlevel > 0)
    options.minimizer_progress_to_stdout = true;
  else
    options.logging_type = ceres::SILENT;
  options.max_linear_solver_iterations = niter;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  if (vlevel > 0)
    std::cout << summary.FullReport() << "\n";

  if (vlevel > 0) {
    int nresid = 3 * (k + attitude.num - 1);
    // total number of residuals
    std::vector<double> resid(nresid);
    ceres::Problem::EvaluateOptions opt;
    opt.apply_loss_function = false;
    double cost;
    problem.Evaluate(opt, &cost, &resid, NULL, NULL);
    double sum = 0;
    for (int i = 0; i < nresid; ++i) sum += resid[i] * resid[i];
    std::cerr << "Initial cost " << cost
              << " RMS resid " << std::sqrt(sum/nresid) << "\n";
  }

  std::chrono::steady_clock::time_point end =
    std::chrono::steady_clock::now();
  std::cerr << "elapsed time "
            << std::chrono::duration_cast<std::chrono::milliseconds>
    (end - start).count()
            << "ms\n";

}

std::pair<Eigen::Vector3d, Eigen::Vector3d> trajfit::coords(double t) const {
  int i = std::min(num-1,
                   std::max(0, int(std::floor((t - tstart) / tblock))));
  double tf = (t - tstart) / tblock - (i + 0.5);
  Eigen::Vector3d pos, vel;
  for (int j = 0; j < 3; ++j)
    pos(j) = lasertraj::EndPointCubic(r[i  ](j), v[i  ](j),
                                      r[i+1](j), v[i+1](j),
                                      tf, vel.data()+j);
  return std::make_pair(pos, vel / tblock);
}

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <GeographicLib/Math.hpp>
#include <LidarTrajectory/PulseCollection.hpp>
#include <ceres/ceres.h>
#include <glog/logging.h>

namespace LidarTrajectory {

  void PulseCollection::SetDefaults() {
    s = START;
    // Require that parameters be set in the configuration file (this makes
    // surprises less likely).
    //
    // dtr = 0.001;
    // dts = 0.001;
    // tblock = 1;
    // dr = 0.01;
    // minsep = 0.1;
    // accelweight = 0.125;
    // clampweight = 0.001;
    // straddleweight = 0.1;
    // vlevel = 0;
    // estn = 20;
    // niter = 50;
    // dang = 1
    // flipscanang = f
    // fixedpitch = nan
    // pitchweight = 1
  }

#define SET(type,val) val = params.lookup<type>(#val)

  void PulseCollection::SetParams() {
    SET(double, dtr);
    SET(double, dts);
    SET(double, tblock);
    SET(double, dr);
    SET(double, minsep);
    SET(double, accelweight);
    SET(double, clampweight);
    SET(double, straddleweight);
    SET(int, vlevel);
    SET(int, estn);
    SET(int, niter);
    SET(bool, flipscanang);
    SET(double, dang);
    SET(double, pitchweight);
    SET(double, scanweightest);
    SET(double, fixedpitch);
    SET(double, multiweight);
    SET(double, scanweight);
    SET(double, attaccelweight);
    SET(double, attclampweight);
    SET(double, extrapitchclamp);
  }

#undef SET

  void PulseCollection::ReadParams(std::istream& str) {
    params.Init(str);
  }
  void PulseCollection::ReadParams(const std::string& file) {
    params.Init(file);
  }
  void PulseCollection::AddParams(std::istream& str) {
    params.Add(str);
  }
  void PulseCollection::AddParams(const std::string& line) {
    params.Add(line);
  }

  void PulseCollection::Add(double t, const Eigen::Vector3d& r,
                            int num, int ret, double ang) {
    // Let's skip return unless num and ret are sensible
    if (!(ret >= 1 && ret <= num))
      return;
    if (flipscanang) ang = -ang;
    if (s == START) {
      Init();
      // Somewhat arbitrary choice of origin (to ensure internal calculations
      // are well conditioned.
      torg = std::floor(t / 10) * 10;
      rorg = (Eigen::floor(r.array() / 100) * 100).matrix();
      s = READING;
      numpending = -1;
      tmin = std::numeric_limits<double>::max(); tmax = -tmin;
    } else if (s == COMPLETE)
      throw std::runtime_error
        ("PulseCollection: Cannot ingest in COMPLETE state");
    /* skipping these for now -- See first line of this function
    if (!(num > 0))
      throw std::runtime_error
        ("PulseCollection: number of returns not positive");
    if (!(ret >= 1 && ret <= num))
      throw std::runtime_error("PulseCollection: invalid return number");
    */
    tmax = std::max(t, tmax);
    tmin = std::min(t, tmin);
    t -= torg;
    // std::cerr << "Add " << std::fixed << std::setprecision(6)
    //           << t << " " << std::setprecision(2)
    //           << r(0) << " "<< r(1) << " " << r(2) << "\n";
    if (numpending < 0) {
      tpending = t;
      rfirst = rlast = r;
      numpending = num;
      angpending = ang;
      minpending = maxpending = ret;
    } else {
      if (t < tpending)
        throw std::runtime_error
          ("PulseCollection: returns are not sorted in time");
      if (t == tpending) {
        if (!(num == numpending && ang == angpending )) {
          // Discard inconsistent returns
          return;
          /*
          std::cerr << std::fixed << std::setprecision(6) << t+torg << " "
                    << num << " " << numpending << " "
                    << ang << " " << angpending << "\n";
          throw std::runtime_error
            ("PulseCollection: inconsistent coincident returns");
          */
        }
        if (ret < minpending) {
          minpending = ret;
          rfirst = r;
        }
        if (ret > maxpending) {
          maxpending = ret;
          rlast = r;
        }
      } else {
        if (std::isfinite(dtr))
          bufferr.push_back(Pulse(tpending,
                                  rfirst - rorg, rlast - rorg,
                                  angpending));
        if (std::isfinite(dts))
          // Mark scan angle pulses with no separation via rlast -> rfirst
          buffers.push_back(Pulse(tpending,
                                  rfirst - rorg, rfirst - rorg,
                                  angpending));
        // std::cerr << "Buffer " << std::fixed << std::setprecision(6)
        //           << tpending << " " << std::setprecision(2)
        //           << rfirst(0) << " "<< rlast(0) << "\n";
        tpending = t;
        rfirst = rlast = r;
        numpending = num;
        angpending = ang;
        minpending = maxpending = ret;
      }
    }
    if (!bufferr.empty() &&
        std::floor(bufferr[0].t / dtr) != std::floor(t / dtr) ) {
      // std::cerr << "Register " << std::fixed << std::setprecision(6)
      //           << buffer[0].t << " " << t << "\n";

      Register(bufferr, true);
    }
    if (!buffers.empty() &&
        std::floor(buffers[0].t / dts) != std::floor(t / dts) ) {
      // std::cerr << "Register " << std::fixed << std::setprecision(6)
      //           << buffer[0].t << " " << t << "\n";

      Register(buffers, false);
    }
  }

  void PulseCollection::Register(std::vector<Pulse>& buffer,
                                 bool multireturn) {
    int num = int(buffer.size());
    if (num == 0) return;
    if (multireturn) {
      int k = -1;
      double d = 0;
      for (int i = 0; i < num; ++i) {
        if (buffer[i].d > d) {
          d = buffer[i].d;
          k = i;
        }
      }
      // Only include pulses if separation > minsep
      if (2*d > minsep) pulses.push_back(buffer[k]);
    } else {
      // Look for midpoint of a run of pulses with the same angle (on the
      // theory that this will have the smallest quantization error).
      int k = num/2,  ka = k - 1, kb = k + 1;
      double ang = buffer[k].ang;
      while (ka >= 0 && buffer[ka].ang == ang) --ka;
      while (kb < num && buffer[kb].ang == ang) ++kb;
      k = (ka + kb) / 2;
      pulses.push_back(buffer[k]);
      // If change in scan angle consistent
      if (ka >= 0 && kb < num &&
          (buffer[k].ang - buffer[ka].ang) *
          (buffer[kb].ang - buffer[k].ang) > 0)
        pulses.back().n = (buffer[k].ang - buffer[ka].ang > 0 ? 1 : -1) *
          (buffer[kb].r - buffer[ka].r);
    }
    buffer.clear();
  }

  void PulseCollection::Finalize() {
    if (s == READING) {
      if (numpending >= 0) {
        if (std::isfinite(dtr)) {
          bufferr.push_back(Pulse(tpending,
                                  rfirst - rorg, rlast - rorg,
                                  angpending));
          Register(bufferr, true);
        }
        if (std::isfinite(dts)) {
          buffers.push_back(Pulse(tpending,
                                  rfirst - rorg, rfirst - rorg,
                                  angpending));
          Register(buffers, false);
        }
        numpending = -1;
      }

      { std::vector<Pulse> z; bufferr.swap(z); }
      { std::vector<Pulse> z; buffers.swap(z); }
    }
    // Don't need to sort the pulses
    // std::sort(pulses.begin(), pulses.end());
    std::string dumppulses = params.lookup<std::string>("dumppulses","");
    if (!dumppulses.empty()) {
      std::ofstream out(dumppulses.c_str());
      out << std::fixed;
      for (const Pulse& p: pulses)
        out << std::setprecision(6)
            << p.t+torg << " "
            << std::setprecision(3)
            << p.r(0)+rorg(0) << " "
            << p.r(1)+rorg(1) << " "
            << p.r(2)+rorg(2) << " "
            << p.d << " "
            << std::setprecision(6)
            << p.n(0) << " " << p.n(1) << " " << p.n(2) << "\n";
    }
    s = COMPLETE;

    Solve();
  }

  bool PulseCollection::EstimatedPositionVelocity(double t,
                                                  Eigen::Vector3d& r,
                                                  Eigen::Vector3d& v)
    const {

    if(vlevel > 1){
      std::cerr<<"EstimatedPositionVelocity for t = "<<std::to_string(t)<<std::endl;
    }

    std::vector<Pulse> psub;
    for (const Pulse& p: pulses) {
      if (p.t >= t - tblock && p.t <= t + tblock) {
        psub.push_back(p);
        psub.back().t -= t;
      }
    }
    int k = int(psub.size());
    if(vlevel > 1){
      std::cerr<<"psub.size() = "<<std::to_string(k)<<std::endl;
    }

    // Accumulate scan vector for scan-angle pulses
    int kscan = 0;
    double sx, sy;
    bool skipscan = false;
    {
      Eigen::Vector3d scandir(Eigen::Vector3d::Zero());
      for (int l = 0; l < k; ++l) {
        if (!psub[l].MultiReturn()) {
          scandir += psub[l].n;
          ++kscan;
        }
      }
      sx = -scandir(0); sy = -scandir(1);
      double h = std::hypot(sx, sy);
      if (h > 0)
        { sx /= h; sy /= h; }
      else
        skipscan = true;
    }

    int m = skipscan ? k - kscan : k;
    if (m  < estn) return false;

    if(vlevel > 1){
      std::cerr<<"rough estimation using m pulses. m="<<std::to_string(m)<<std::endl;
    }

    // For multi-return pulses
    // equation for pulse starting at r in direction n, distance = s
    //
    //   x + vx*t = rx + nx * s
    //   y + vy*t = ry + ny * s
    //   z + vz*t = rz + nz * s
    //
    // replace s by z as parameterization, s = (z + t*vz - rz)/nz
    //   nz*x      - nx*z + nz*t*vx           - nx*t*vz = nz*rx - nx*rz
    //        nz*y - ny*z           + nz*t*vy - ny*t*vz = nz*ry - ny*rz
    // or
    //
    //    A . [x,y,z,vx,vy,vz]' = B
    //
    // where
    //
    //   A = [nz,  0, -nx, nz*t,    0, -nx*t]
    //       [ 0, nz, -ny,    0, nz*t, -ny*t]
    //       ... for all pulses
    //   B = [nz * rx - nx * rz]
    //       [nz * ry - ny * rz]
    //       ... for all pulses
    // finally, we'll weight use pulse contribution by d

    // For scan angle pulses, we can substitute
    //
    //   nx = sin(ang) * sx, ny = sin(ang) * sy, nz = cos(ang)
    //
    // Thus assumes pitch = 0.  It would be simple to use fixedpitch here.
    Eigen::MatrixXd A(2 * m, 6);
    Eigen::VectorXd B(2 * m);
    for (int l = 0, h = 0; l < k; ++l) {
      const Pulse& p = psub[l];
      double nx, ny, nz, w;
      if (p.MultiReturn()) {
        nx = p.n(0); ny = p.n(1); nz = p.n(2);
        w = p.d;
      } else {
        if (skipscan) continue;
        double sang = sin(p.ang), cang = cos(p.ang);
        nx = sang * sx; ny = sang * sy; nz = cang;
        w = scanweightest;
      }
      A(2*h+0, 0) = nz; A(2*h+0, 2) = -nx;
      A(2*h+1, 1) = nz; A(2*h+1, 2) = -ny;
      A(2*h+0, 3) = nz*p.t; A(2*h+0, 5) = -nx*p.t;
      A(2*h+1, 4) = nz*p.t; A(2*h+1, 5) = -ny*p.t;
      A(2*h+0, 1) = A(2*h+0, 4) = A(2*h+1, 0) = A(2*h+1, 3) = 0;
      B(2*h+0) = nz * p.r(0) - nx * p.r(2);
      B(2*h+1) = nz * p.r(1) - ny * p.r(2);
      A.row(2*h+0) *= w; B(2*h+0) *= w;
      A.row(2*h+1) *= w; B(2*h+1) *= w;
      ++h;
    }
    Eigen::Matrix<double, 1, Eigen::Dynamic> rv(A.jacobiSvd(Eigen::ComputeThinU |
                                               Eigen::ComputeThinV).solve(B));
    r = rv.head<3>();
    v = rv.tail<3>();
    return true;
  }

  void PulseCollection::InitializeTrajectory() {

    if(vlevel>0)
    {
      std::cerr<<std::endl<<"InitializeTrajectory() "<<std::endl;
    }

    int npulses = int(pulses.size());

    if(vlevel>0)
    {
      std::cerr<<"num pulses: "<<std::to_string(npulses)<<std::endl;
    }

    if (npulses == 0)
      throw std::runtime_error("PulseCollection: no pulses for Solve");
    double tstart = std::floor((tmin - torg) / tblock);
    int num = int(std::ceil((tmax - torg) / tblock) - tstart) - 1;
    if (num < 1)
      throw std::runtime_error("PulseCollection: no time interval for Solve");
    tstart *= tblock;

    if(vlevel>0)
    {
      std::cerr<<"num traj/attitude: "<<std::to_string(num)<<std::endl;
    }
    traj = SplineFit3(num, tblock, tstart);
    attitude = SplineFit2(num, tblock, tstart);

    for (int i = 0; i <= num; ++i) {
      traj.missing[i] = !EstimatedPositionVelocity(tstart + i * tblock,
                                                   traj.r[i], traj.v[i]);
      traj.v[i] *= tblock;
    }

    if (!traj.fillmissing(true))
      throw std::runtime_error
        ("PulseCollection: to few pulses for initial estimate of trajectory");

    for (int i = 0; i <= num; ++i) {
      int im = std::max(0, i-1), ip = std::min(num, i+1);
      // atan(dx, dy) to give clockwise from north convention
      attitude.r[i] =
        Eigen::Vector2d(std::atan2(traj.r[ip](0) - traj.r[im](0),
                                   traj.r[ip](1) - traj.r[im](1)),
                        std::isnan(fixedpitch) ? 0.0 :
                        fixedpitch * GeographicLib::Math::degree());
      attitude.v[i] = Eigen::Vector2d::Zero();
    }
    // Make sure heading doesn't jump around
    double ang0 = attitude.r[0](0) / GeographicLib::Math::degree();
    for (int i = 1; i <= num; ++i) {
      double ang1 = attitude.r[i](0) / GeographicLib::Math::degree();
      // AngDiff(x, y) returns y - x in [-180, 180]
      attitude.r[i](0) = (ang0 + GeographicLib::Math::AngDiff(ang0, ang1)) *
        GeographicLib::Math::degree();
    }

    if(vlevel>0)
    {
      std::cerr<<"completed traj/attitude initialization: EstimatedPositionVelocity()."<<std::endl;
    }

    std::string dumpinittraj = params.lookup<std::string>("dumpinittraj","");
    if (!dumpinittraj.empty()) {
      std::ofstream str(dumpinittraj.c_str());
      str << std::fixed << std::setprecision(3);
      for (int i = 0; i <= num; ++i)
        str << tstart + i * tblock + torg << " "
            << traj.r[i](0) + rorg(0) << " "
            << traj.r[i](1) + rorg(1) << " "
            << traj.r[i](2) + rorg(2) << " "
            << traj.v[i](0) / tblock << " "
            << traj.v[i](1) / tblock << " "
            << traj.v[i](2) / tblock << " "
            << attitude.r[i](0) / GeographicLib::Math::degree() << " "
            << attitude.r[i](1) / GeographicLib::Math::degree() << " "
            << attitude.v[i](0) / GeographicLib::Math::degree() / tblock << " "
            << attitude.v[i](1) / GeographicLib::Math::degree() / tblock << " "
            << traj.missing[i] << "\n";
    }
    if(vlevel>0)
    {
      std::cerr<<"completed initializing traj and attitude "<<std::to_string(npulses)<<std::endl;
    }
  }

  void PulseCollection::Solve() {
    if(vlevel>0)
    {
      std::cerr<<std::endl<<"Solve()"<<std::endl;
    }

    InitializeTrajectory();
    int num = traj.num;
    if(vlevel>0)
    {
      std::cerr<<"num traj: "<<std::to_string(num)<<std::endl;
    }

    google::InitGoogleLogging("LidarTrajectory");
    ceres::Problem problem;

    int npulses = int(pulses.size());
    if(vlevel>0)
    {
      std::cerr<<"num pulses: "<<std::to_string(npulses)<<std::endl;
      std::cerr<<"num traj: "<<std::to_string(num)<<std::endl;
    }

    // for debugging
    int numMultiReturnResidualBlocks = 0;
    int numScanAngleResidualBlocks = 0;

    // Set up residual blocks for pulses
    for (int j = 0; j < npulses; ++j) {
      const Pulse& p = pulses[j];
      auto tconv = traj.tconvert(p.t);
      int i = tconv.first;
      double t = tconv.second;
      if (p.MultiReturn()) {
        ceres::CostFunction* cost_function =
          new ceres::AutoDiffCostFunction<FirstLastError,
                                          2,       // number of residuals
                                          3,3,3,3> // data for cubic fit
          (new FirstLastError(p, t));
        ceres::LossFunction* loss_function =
          new ceres::ScaledLoss(new ceres::CauchyLoss(dr),
                                multiweight,
                                ceres::TAKE_OWNERSHIP);
        problem.AddResidualBlock(cost_function, loss_function,
                                 traj.r[i  ].data(),
                                 traj.v[i  ].data(),
                                 traj.r[i+1].data(),
                                 traj.v[i+1].data());
        numMultiReturnResidualBlocks++;
      } else {
        ceres::CostFunction* cost_function =
          new ceres::AutoDiffCostFunction<ScanAngleError,
                                          2,       // number of residuals
                                          3,3,3,3,
                                          2,2,2,2> // data for cubic fit
          (new ScanAngleError(p, t, pitchweight,
                              fixedpitch * GeographicLib::Math::degree()));
        ceres::LossFunction* loss_function = new ceres::ScaledLoss
          (new ceres::CauchyLoss(dang * GeographicLib::Math::degree()),
           scanweight,
           ceres::TAKE_OWNERSHIP);
        problem.AddResidualBlock(cost_function, loss_function,
                                 traj.r[i  ].data(),
                                 traj.v[i  ].data(),
                                 traj.r[i+1].data(),
                                 traj.v[i+1].data(),
                                 attitude.r[i  ].data(),
                                 attitude.v[i  ].data(),
                                 attitude.r[i+1].data(),
                                 attitude.v[i+1].data());
        numScanAngleResidualBlocks++;
      }
    }

    if(vlevel>0)
    {
      std::cerr<<"num multi-return residual blocks: "<<std::to_string(numMultiReturnResidualBlocks)<<std::endl;
      std::cerr<<"num multi-return residual blocks: "<<std::to_string(numScanAngleResidualBlocks)<<std::endl;
    }

    // The acceleration constraints for traj and attitude
    if (accelweight > 0) {
      for (int i = 1; i < num; ++i) {
        ceres::CostFunction* cost_function =
          new ceres::AutoDiffCostFunction<AccelJumpConstraint<3>,
                                          3,       // number of residuals
                                          3,3,3,3,3> // data for cubic fit
          (new AccelJumpConstraint<3>(tblock));
        ceres::LossFunction* loss_function =
          (ceres::LossFunction*)
          (new ceres::ScaledLoss(nullptr, accelweight,
                                 ceres::TAKE_OWNERSHIP));
        problem.AddResidualBlock(cost_function, loss_function,
                                 traj.r[i-1].data(),
                                 traj.v[i-1].data(),
                                 traj.v[i  ].data(),
                                 traj.r[i+1].data(),
                                 traj.v[i+1].data());
      }
    }
    if (attaccelweight > 0) {
      for (int i = 1; i < num; ++i) {
        ceres::CostFunction* cost_function =
          new ceres::AutoDiffCostFunction<AccelJumpConstraint<2>,
                                          2,       // number of residuals
                                          2,2,2,2,2> // data for cubic fit
          (new AccelJumpConstraint<2>(tblock));
        ceres::LossFunction* loss_function =
          (ceres::LossFunction*)
          (new ceres::ScaledLoss(nullptr, attaccelweight,
                                 ceres::TAKE_OWNERSHIP));
        problem.AddResidualBlock(cost_function, loss_function,
                                 attitude.r[i-1].data(),
                                 attitude.v[i-1].data(),
                                 attitude.v[i  ].data(),
                                 attitude.r[i+1].data(),
                                 attitude.v[i+1].data());
      }
    }

    // Also allow "clamping" of at the nodes with constrains neighboring cubic
    // polynomials to be close to one another.  In general clampweight should
    // be "small".  straddleweight is a larger weight to enforce clamping where
    // there's a sparsity of data.
    if (clampweight > 0 || straddleweight > 0) {
      for (int i = 1; i < num; ++i) {
        double w = traj.missing[i] ? straddleweight : clampweight;
        if (!(w > 0)) continue;
        ceres::CostFunction* cost_function =
          new ceres::AutoDiffCostFunction<ClampConstraint<3>,
                                          3,       // number of residuals
                                          3,3,3,3,3> // data for cubic fit
          (new ClampConstraint<3>(tblock));
        ceres::LossFunction* loss_function =
          (ceres::LossFunction*)
          (new ceres::ScaledLoss(nullptr, w,
                                 ceres::TAKE_OWNERSHIP));
        problem.AddResidualBlock(cost_function, loss_function,
                                 traj.r[i-1].data(),
                                 traj.v[i-1].data(),
                                 traj.r[i  ].data(),
                                 traj.r[i+1].data(),
                                 traj.v[i+1].data());
      }
    }
    // The estimate of the pitch can sometimes oscillate too much.
    // extrapitchclamp is a way to suppress this.
    if (attclampweight > 0) {
      Eigen::Vector2d mult(1.0, extrapitchclamp);
      for (int i = 1; i < num; ++i) {
        ceres::CostFunction* cost_function =
          new ceres::AutoDiffCostFunction<ClampConstraint<2>,
                                          2,       // number of residuals
                                          2,2,2,2,2> // data for cubic fit
          (new ClampConstraint<2>(tblock, mult));
        ceres::LossFunction* loss_function =
          (ceres::LossFunction*)
          (new ceres::ScaledLoss(nullptr, attclampweight,
                                 ceres::TAKE_OWNERSHIP));
        problem.AddResidualBlock(cost_function, loss_function,
                                 attitude.r[i-1].data(),
                                 attitude.v[i-1].data(),
                                 attitude.r[i  ].data(),
                                 attitude.r[i+1].data(),
                                 attitude.v[i+1].data());
      }
    }

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
      std::cerr << summary.FullReport() << "\n";

    // for debug - print out few lines of traj and attitude
    if(vlevel>0)
    {
      std::cerr<<"Estimated traj and attitude (first few values)"<<std::endl;
      for (int i = 0; i < std::min(10, traj.num); ++i) {
        std::cerr<<" traj "<<std::to_string(i)<<" ("
                           <<std::to_string(traj.r[i](0))<<", "
                           <<std::to_string(traj.r[i](1))<<", "
                           <<std::to_string(traj.r[i](2))<<") ";
        std::cerr<<" attitude "<<std::to_string(i)<<" ("
                               <<std::to_string(attitude.r[i](0))<<", "
                               <<std::to_string(attitude.r[i](1))<<")"<<std::endl;
      }
    }
  }

  Eigen::Vector3d
  PulseCollection::Trajectory(double t, Eigen::Vector3d& v, Eigen::Vector3d& a)
    const {
    if (s != COMPLETE)
      throw std::runtime_error("PulseCollection: not yet finalized");
    Eigen::Vector3d pos = traj.position(t - torg, v, a) + rorg;
    return pos;
  }

  Eigen::Vector2d
  PulseCollection::Attitude(double t, Eigen::Vector2d& v)
    const {
    if (s != COMPLETE)
      throw std::runtime_error("PulseCollection: not yet finalized");
    Eigen::Vector2d p = attitude.position(t - torg, v);
    p /= GeographicLib::Math::degree();
    v /= GeographicLib::Math::degree();
    return p;
  }

  std::string PulseCollection::ConfigRoot() {
    char* configroot = getenv("SRI_TRAJECTORY_CONFIG_DIR");
    std::string DIR = configroot ? std::string(configroot) : ".";
    return DIR;
  }

}

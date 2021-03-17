#if !defined(LIDARTRAJECTORY_PULSECOLLECTION_HPP)
#define LIDARTRAJECTORY_PULSECOLLECTION_HPP 1

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <LidarTrajectory/Pulse.hpp>
#include <LidarTrajectory/SplineFit.hpp>
#include <LidarTrajectory/Params.hpp>

namespace LidarTrajectory {

  class PulseCollection {
    // Samples and saves lidar pulses
  private:
    enum state { START, READING, COMPLETE };
    state s;

    std::vector<Pulse> bufferr, buffers;
    double tpending, angpending;
    int numpending, minpending, maxpending;
    Eigen::Vector3d rfirst, rlast;
    // Select a pulse to save
    void Register(std::vector<Pulse>& buffer, bool multireturn);
    // Estimate pos + vel at t based on pulses within +/- tblock at t.
    bool EstimatedPositionVelocity(double t, Eigen::Vector3d& r,
                                   Eigen::Vector3d& v) const;
    void InitializeTrajectory();
    // Do the ceres solution
    void Solve();
    void SetDefaults();
    void SetParams();
  public:
    enum mode { MULTIRETURNS = 0, ATTITUDE, ATTITUDEFULL };
    Params params;

    double torg, dtr, dts, tblock, dr, minsep, multiweight,
      accelweight, clampweight, straddleweight, tmin, tmax;
    int vlevel,estn, niter;
    Eigen::Vector3d rorg;
    std::vector<Pulse> pulses;
    SplineFit3 traj;

    // scan angle
    double dang, fixedpitch, pitchweight, scanweightest, scanweight,
      attaccelweight, attclampweight, extrapitchclamp;
    bool flipscanang;
    SplineFit2 attitude;

    void Init() {
      SetParams();
      torg = 0;
      rorg = Eigen::Vector3d::Zero();
      { std::vector<Pulse> z; bufferr.swap(z); }
      { std::vector<Pulse> z; buffers.swap(z); }
      { std::vector<Pulse> z; pulses.swap(z); }
      s = START;
      numpending = -1;          // Use as flag when no there's no pending data
    }
    PulseCollection() {
      SetDefaults();
      std::string config = ConfigRoot() + "/trajectory.cfg";
      std::ifstream f(config.c_str());
      if (f.good())
        ReadParams(f);
      else
        throw std::runtime_error("Cannot open configuration file " + config);
    }
    PulseCollection(std::istream& str) {
      SetDefaults(); ReadParams(str);
    }
    PulseCollection(const std::string& file) {
      SetDefaults(); ReadParams(file);
    }
    void ReadParams(std::istream& str);
    void ReadParams(const std::string& file);
    void AddParams(std::istream& str);
    void AddParams(const std::string& line);
    // Ingest a lidar return
    void Add(double t, const Eigen::Vector3d& r, int num, int ret, double ang);
    // Finalize a collection and find the trajectory
    void Finalize();
    Eigen::Vector3d Trajectory(double t,
                               Eigen::Vector3d& v, Eigen::Vector3d& a) const;
    Eigen::Vector2d Attitude(double t, Eigen::Vector2d& v) const;
    /* Return the environment variable SRI_QTM_CONFIG_DIR.  If it's
     * not defined return ".". */
    static std::string ConfigRoot();
  };

}
#endif

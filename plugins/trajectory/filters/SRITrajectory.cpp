// SRITrajectory.cpp

#include "SRITrajectory.hpp"
#include <iostream>
#include <limits>
#include <cmath>
#include <Eigen/Dense>

#include <LidarTrajectory/PulseCollection.hpp>

namespace pdal {

  using namespace std;

  static PluginInfo const s_info
    {
     "filters.sritrajectory",
     "SRI Trajectory calculator",
     "http://link/to/documentation"
    };

struct PrivateArgs {
  StringList m_params;
};

  CREATE_SHARED_STAGE(SRITrajectory, s_info)

  SRITrajectory::SRITrajectory()
    : Filter(),
    m_args(new PrivateArgs)
  {}

  string SRITrajectory::getName() const { return s_info.name; }

  void SRITrajectory::addArgs(ProgramArgs& args) {
    args.add("params", "Set params on backend engine using "
             "\"key1 = val1, key2 = val2, ...\".  Useful parameters "
             "(and defaults) are: dt (0.0005 s) the sampling interval; "
             "tblock (1 s) the spline block size; "
             "tout (0.005 s) the time resolustion for the output trajectory; "
             "dr (0.01 m) the resolution of the lidar positions.",
             m_args->m_params, StringList(0));
  }

  void SRITrajectory::addDimensions(PointLayoutPtr layout) {
    layout->registerDim(Dimension::Id::X);
    layout->registerDim(Dimension::Id::Y);
    layout->registerDim(Dimension::Id::Z);
    layout->registerDim(Dimension::Id::XVelocity);
    layout->registerDim(Dimension::Id::YVelocity);
    layout->registerDim(Dimension::Id::ZVelocity);
    layout->registerDim(Dimension::Id::XBodyAccel);
    layout->registerDim(Dimension::Id::YBodyAccel);
    layout->registerDim(Dimension::Id::ZBodyAccel);
    layout->registerDim(Dimension::Id::Azimuth);
    // layout->registerDim(Dimension::Id::Roll); not set
    layout->registerDim(Dimension::Id::Pitch);
    layout->registerDim(Dimension::Id::XBodyAngRate); // dPitch/dt
    // layout->registerDim(Dimension::Id::YBodyAngRate); not set
    layout->registerDim(Dimension::Id::ZBodyAngRate); // dAzimuth/dt
  }

  void SRITrajectory::prepared(PointTableRef /* table */)  {}

  PointViewSet SRITrajectory::run(PointViewPtr inView) {
    using namespace pdal;
    PointViewPtr outView = inView->makeNew();
    // Segment input view into ignored/kept views.
    PointViewPtr ignoredView = inView->makeNew();
    PointViewPtr keptView = inView->makeNew();
    keptView->append(*inView);
    PointViewPtr syntheticView = keptView->makeNew();
    PointViewPtr realView = keptView->makeNew();
    realView->append(*keptView);

    // Segment kept view into two views
    PointViewPtr inlierView = realView->makeNew();
    inlierView->append(*realView);

    if (!inlierView->size())
    {
        throwError("No returns to process.");
    }

    LidarTrajectory::PulseCollection coll;
    for (const string& p: m_args->m_params)
      coll.params.Add(p);

    log()->get(LogLevel::Debug2) << "flipscanang = "<<coll.flipscanang<<std::endl;

    Eigen::Vector3d r;
    for (PointId i = 0; i < inlierView->size(); ++i) {
      PointRef point = inlierView->point(i);
      double t = point.getFieldAs<double>(Dimension::Id::GpsTime);
      r(0) = point.getFieldAs<double>(Dimension::Id::X);
      r(1) = point.getFieldAs<double>(Dimension::Id::Y);
      r(2) = point.getFieldAs<double>(Dimension::Id::Z);
      double ang = point.getFieldAs<double>(Dimension::Id::ScanAngleRank);
      int num = int(point.getFieldAs<unsigned char>
                    (Dimension::Id::NumberOfReturns));
      int ret = int(point.getFieldAs<unsigned char>
                    (Dimension::Id::ReturnNumber));
      coll.Add(t, r, num, ret, ang);
    }

    coll.Finalize();

    double tout = coll.params.lookup<double>("tout", 0.005);
    double tmin = floor(coll.tmin / tout);
    int nt = int(ceil(coll.tmax / tout) - tmin);
    tmin *= tout;

    Eigen::Vector3d pos, vel, accel;
    Eigen::Vector2d att, attvel;
    for (int it = 0; it <= nt; ++it) {
      double t = tmin + it * tout;
      pos = coll.Trajectory(t, vel, accel);
      att = coll.Attitude(t, attvel);

      PointId idx = outView->size();
      outView->setField(Dimension::Id::GpsTime,      idx, t);
      outView->setField(Dimension::Id::X,            idx, pos(0));
      outView->setField(Dimension::Id::Y,            idx, pos(1));
      outView->setField(Dimension::Id::Z,            idx, pos(2));
      outView->setField(Dimension::Id::XVelocity,    idx, vel(0));
      outView->setField(Dimension::Id::YVelocity,    idx, vel(1));
      outView->setField(Dimension::Id::ZVelocity,    idx, vel(2));
      outView->setField(Dimension::Id::XBodyAccel,   idx, accel(0));
      outView->setField(Dimension::Id::YBodyAccel,   idx, accel(1));
      outView->setField(Dimension::Id::ZBodyAccel,   idx, accel(2));
      outView->setField(Dimension::Id::Azimuth,      idx, att(0));
      outView->setField(Dimension::Id::Pitch,        idx, att(1));
      outView->setField(Dimension::Id::ZBodyAngRate, idx, attvel(0));
      outView->setField(Dimension::Id::XBodyAngRate, idx, attvel(1));
    }

    log()->get(LogLevel::Debug1) << "Completed trajectory\n";

    PointViewSet viewSet;
    viewSet.insert(outView);
    return viewSet;
  }

} // namespace pdal

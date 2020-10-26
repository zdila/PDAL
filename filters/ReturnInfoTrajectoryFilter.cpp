/******************************************************************************
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of Hobu, Inc. or Flaxen Geo Consulting nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 ****************************************************************************/

#include "ReturnInfoTrajectoryFilter.hpp"

#include "private/trajectory/PulseCollection.hpp"

#include <Eigen/Dense>

namespace pdal
{

using namespace Dimension;
using namespace Eigen;

static StaticPluginInfo const s_info{
    "filters.returninfotrajectory", "Return Info Trajectory estimation",
    "htts://pdal.io/stages/filters.returninfotrajectory.html"};

struct PrivateArgs
{
    double m_dt;
    double m_tblock;
    double m_dr;
    double m_minsep;
    double m_accelweight;
    double m_clampweight;
    double m_straddleweight;
    int m_vlevel;
    int m_estn;
    int m_niter;
    double m_tout;
};

CREATE_STATIC_STAGE(ReturnInfoTrajectory, s_info)

ReturnInfoTrajectory::ReturnInfoTrajectory() : Filter(), m_args(new PrivateArgs)
{
}

std::string ReturnInfoTrajectory::getName() const
{
    return s_info.name;
}

void ReturnInfoTrajectory::addArgs(ProgramArgs& args)
{
    args.add("dt", "Sampling interval (s)", m_args->m_dt, 0.001);
    args.add("tblock", "Block size for cubic spline (s)", m_args->m_tblock,
             1.0);
    args.add("dr", "Error in returns (m)", m_args->m_dr, 0.01);
    args.add("minsep", "Min separation of returns considered", m_args->m_minsep,
             0.1);
    args.add("accelweight", "Weighting for acceleration constraint",
             m_args->m_accelweight, 0.125);
    args.add("clampweight", "Weighting for general clamping constraint",
             m_args->m_clampweight, 0.001);
    args.add("straddleweight", "Weighting for straddling clamping constraint",
             m_args->m_straddleweight, 0.1);
    args.add("vlevel", "Verbosity level", m_args->m_vlevel, 0);
    args.add("estn", "Min number of pulses to use of rough estimate",
             m_args->m_estn, 20);
    args.add("niter", "Number of iterations for Ceres solver", m_args->m_niter,
             50);
    args.add("tout", "Interval for the reported trajectory", m_args->m_tout,
             0.01);
}

void ReturnInfoTrajectory::addDimensions(PointLayoutPtr layout)
{
    layout->registerDim(Id::X);
    layout->registerDim(Id::Y);
    layout->registerDim(Id::Z);
    layout->registerDim(Id::XVelocity);
    layout->registerDim(Id::YVelocity);
    layout->registerDim(Id::ZVelocity);
    layout->registerDim(Id::XBodyAccel);
    layout->registerDim(Id::YBodyAccel);
    layout->registerDim(Id::ZBodyAccel);
}

void ReturnInfoTrajectory::prepared(PointTableRef table) {}

PointViewSet ReturnInfoTrajectory::run(PointViewPtr inView)
{
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

    LidarTrajectory::PulseCollection coll(
        m_args->m_dt, m_args->m_tblock, m_args->m_dr, m_args->m_minsep,
        m_args->m_accelweight, m_args->m_clampweight, m_args->m_straddleweight,
        m_args->m_vlevel, m_args->m_estn, m_args->m_niter);
    Vector3d r;
    for (PointRef point : *inlierView)
    {
        double t = point.getFieldAs<double>(Id::GpsTime);
        r(0) = point.getFieldAs<double>(Id::X);
        r(1) = point.getFieldAs<double>(Id::Y);
        r(2) = point.getFieldAs<double>(Id::Z);
        double ang = point.getFieldAs<double>(Id::ScanAngleRank);
        int num = int(point.getFieldAs<unsigned char>(Id::NumberOfReturns));
        int ret = int(point.getFieldAs<unsigned char>(Id::ReturnNumber));
        coll.Add(t, r, num, ret, ang);
    }
    coll.Finalize();
    double tmin = floor(coll.tmin / m_args->m_tout);
    int nt = int(ceil(coll.tmax / m_args->m_tout) - tmin);
    tmin *= m_args->m_tout;

    Vector3d pos, vel, accel;
    for (int it = 0; it <= nt; ++it)
    {
        double t = tmin + it * m_args->m_tout;
        pos = coll.Trajectory(t, vel, accel);
        PointId idx = outView->size();
        outView->setField(Id::GpsTime, idx, t);
        outView->setField(Id::X, idx, pos(0));
        outView->setField(Id::Y, idx, pos(1));
        outView->setField(Id::Z, idx, pos(2));
        outView->setField(Id::XVelocity, idx, vel(0));
        outView->setField(Id::YVelocity, idx, vel(1));
        outView->setField(Id::ZVelocity, idx, vel(2));
        outView->setField(Id::XBodyAccel, idx, accel(0));
        outView->setField(Id::YBodyAccel, idx, accel(1));
        outView->setField(Id::ZBodyAccel, idx, accel(2));
    }

    log()->get(LogLevel::Debug1) << "Completed trajectory\n";

    PointViewSet viewSet;
    viewSet.insert(outView);
    return viewSet;
}

} // namespace pdal

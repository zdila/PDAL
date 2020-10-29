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

#pragma once

#include "Pulse.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace LidarTrajectory
{

using namespace Eigen;

class PulseCollection
{
    // Samples and saves lidar pulses
private:
    enum state
    {
        START,
        READING,
        COMPLETE
    };
    state s;
    std::vector<Pulse> buffer;
    double tpending, angpending;
    int numpending, minpending, maxpending;
    Vector3d rfirst, rlast;
    // Select a pulse to save
    void Register();
    // Estimate position at t based on estn multi-return pulses closest to t.
    Vector3d EstimatedPosition(double t, double& estt) const;
    // Estimate pos + vel at t based on pulses within +/- tblock at t.
    bool EstimatedPositionVelocity(double t, Vector3d& r, Vector3d& v) const;
    // Do the ceres solution
    void Solve();
    void SetDefaults();

public:
    double torg, m_dt, m_tblock, m_dr, m_minsep, m_accelweight, m_clampweight,
        m_straddleweight, tmin, tmax;
    int m_vlevel, m_estn, m_niter;
    Vector3d rorg;
    std::vector<Pulse> pulses;
    SplineFit3 traj;
    void Init()
    {
        torg = 0;
        rorg = Vector3d::Zero();
        {
            std::vector<Pulse> z;
            buffer.swap(z);
        }
        {
            std::vector<Pulse> z;
            pulses.swap(z);
        }
        s = START;
        numpending = -1; // Use as flag when no there's no pending data
    }
    PulseCollection(double dt, double tblock, double dr, double minsep,
                    double accelweight, double clampweight,
                    double straddleweight, int vlevel, int estn, int niter)
        : m_dt(dt), m_tblock(tblock), m_dr(dr), m_minsep(minsep),
          m_accelweight(accelweight), m_clampweight(clampweight),
          m_straddleweight(straddleweight), m_vlevel(vlevel), m_estn(estn),
          m_niter(niter)
    {
        SetDefaults();
    }
    // Ingest a lidar return
    void Add(double t, const Vector3d& r, int num, int ret, double ang);
    // Finalize a collection and find the trajectory
    void Finalize();
    Vector3d Trajectory(double t, Vector3d& v, Vector3d& a) const;
};

} // namespace LidarTrajectory

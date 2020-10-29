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

#include "SplineFit.hpp"
#include "Utils.hpp"

#include <Eigen/Dense>

namespace LidarTrajectory
{

using namespace Eigen;

class Pulse
{
    // Data on one pulse of the lidar.  If multiple returns, gives the midpoint
    // of the first/last returns and the unit vector from last to first.  d is
    // the distance from the midpoint to the first return.
public:
    double t;
    Vector3d r, n;
    double d, ang;
    Pulse() : t(0), r(Vector3d::Zero()), n(Vector3d::Zero()), d(0), ang(0) {}
    Pulse(double _t, const Vector3d& _r, double _ang = 0);
    Pulse(double _t, const Vector3d& _f, const Vector3d& _l, double _ang = 0);
    bool operator<(const Pulse& p) const
    {
        return t < p.t;
    }
};

class FirstLastError
{
    double _t;
    const Vector3d _r;
    Matrix3d _M;

public:
    FirstLastError(const Pulse& p, double t)
        : _t(t) // fractional time in block
          ,
          _r(p.r), _M(Utils::PerpProjector(p.n, p.d))
    {
    }
    template <typename T>
    bool operator()(const T* const rm, // 3 vec for pos at beg
                    const T* const vm, // 3 vec for vel at beg
                    const T* const rp, // 3 vec for pos at end
                    const T* const vp, // 3 vec for vel at end
                    // 2 residuals
                    T* residuals) const
    {
        const T t = T(_t), r[3] = {T(_r(0)), T(_r(1)), T(_r(2))};
        T p[3]; // Platform relative to return
        for (int i = 0; i < 3; ++i)
            p[i] =
                SplineFitScalar::EndPointCubic(rm[i], vm[i], rp[i], vp[i], t) -
                r[i];
        // Transform
        T xt = T(_M(0, 0)) * p[0] + T(_M(0, 1)) * p[1] + T(_M(0, 2)) * p[2];
        T yt = T(_M(1, 0)) * p[0] + T(_M(1, 1)) * p[1] + T(_M(1, 2)) * p[2];
        T zt = T(_M(2, 0)) * p[0] + T(_M(2, 1)) * p[1] + T(_M(2, 2)) * p[2];
        // Project
        residuals[0] = xt / zt;
        residuals[1] = yt / zt;
        return true;
    }
};

} // namespace LidarTrajectory

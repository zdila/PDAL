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

#include <Eigen/Dense>

namespace LidarTrajectory
{

using namespace Eigen;

class Utils
{
public:
    static Matrix3d PerpProjector(const Vector3d& v, double d)
    {
        // Assume v is a unit vector
        double v2 = v(0) * v(0) + v(1) * v(1);
        Matrix3d M;
        // See proj.mac for derivation
        M(0, 0) = (v(0) * v(0) * v(2) + v(1) * v(1)) / v2;
        M(1, 1) = (v(1) * v(1) * v(2) + v(0) * v(0)) / v2;
        M(0, 1) = -(1 - v(2)) * v(0) * v(1) / v2;
        M(1, 0) = M(0, 1);
        M(0, 2) = -v(0);
        M(1, 2) = -v(1);
        M.row(2) = v.transpose();
        // Scale first two rows to get projection error at v
        M.row(0) *= d;
        M.row(1) *= d;
        return M;
    }
};

} // namespace LidarTrajectory

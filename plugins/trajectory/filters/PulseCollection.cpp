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

#include "PulseCollection.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <glog/logging.h>

namespace LidarTrajectory
{

void PulseCollection::SetDefaults()
{
    s = START;
    // Require that parameters be set in the configuration file (this makes
    // surprises less likely).
    //
    // m_dt = 0.001;
    // m_tblock = 1;
    // m_dr = 0.01;
    // m_minsep = 0.1;
    // m_accelweight = 0.125;
    // m_clampweight = 0.001;
    // m_straddleweight = 0.1;
    // m_vlevel = 0;
    // m_estn = 20;
    // m_niter = 50;
}

void PulseCollection::Add(double t, const Vector3d& r, int num, int ret,
                          double ang)
{
    // Let's skip return unless num and ret are sensible
    if (!(ret >= 1 && ret <= num))
        return;
    if (s == START)
    {
        Init();
        // Somewhat arbitrary choice of origin (to ensure internal calculations
        // are well conditioned.
        torg = std::floor(t / 10) * 10;
        rorg = (Eigen::floor(r.array() / 100) * 100).matrix();
        s = READING;
        numpending = -1;
        tmin = std::numeric_limits<double>::max();
        tmax = -tmin;
    }
    else if (s == COMPLETE)
        throw std::runtime_error(
            "PulseCollection: Cannot ingest in COMPLETE state");
    /* skipping these for now -- See first line of this function
    if (!(num > 0))
      throw std::runtime_error
        ("PulseCollection: number of returns not positive");
    if (!(ret >= 1 && ret <= num))
      throw std::runtime_error("PulseCollection: invalid return number");
    */
    if (num <= 1)
        return;
    tmax = std::max(t, tmax);
    tmin = std::min(t, tmin);
    t -= torg;
    // std::cerr << "Add " << std::fixed << std::setprecision(6)
    //           << t << " " << std::setprecision(2)
    //           << r(0) << " "<< r(1) << " " << r(2) << "\n";
    if (numpending < 0)
    {
        tpending = t;
        rfirst = rlast = r;
        numpending = num;
        angpending = ang;
        minpending = maxpending = ret;
    }
    else
    {
        if (t < tpending)
            throw std::runtime_error(
                "PulseCollection: returns are not sorted in time");
        if (t == tpending)
        {
            if (!(num == numpending && ang == angpending))
            {
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
            if (ret < minpending)
            {
                minpending = ret;
                rfirst = r;
            }
            if (ret > maxpending)
            {
                maxpending = ret;
                rlast = r;
            }
        }
        else
        {
            buffer.push_back(
                Pulse(tpending, rfirst - rorg, rlast - rorg, angpending));
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
    if (!buffer.empty() &&
        std::floor(buffer[0].t / m_dt) != std::floor(t / m_dt))
    {
        // std::cerr << "Register " << std::fixed << std::setprecision(6)
        //           << buffer[0].t << " " << t << "\n";

        Register();
    }
}

void PulseCollection::Register()
{
    int num = int(buffer.size()), k = -1;
    if (num == 0)
        return;
    double d = -1;
    for (int i = 0; i < num; ++i)
    {
        if (buffer[i].d > d)
        {
            d = buffer[i].d;
            k = i;
        }
    }
    // Use "middle" pulse if no multiple returns;
    // if (d == 0) k = num/2;
    // Only includer pulses if separation > m_minsep
    if (2 * d > m_minsep)
        pulses.push_back(buffer[k]);
    buffer.clear();
}

void PulseCollection::Finalize()
{
    if (s == READING)
    {
        if (numpending >= 0)
        {
            buffer.push_back(
                Pulse(tpending, rfirst - rorg, rlast - rorg, angpending));
            Register();
            numpending = -1;
        }
    }
    // Don't need to sort the pulses
    // std::sort(pulses.begin(), pulses.end());
    std::string dumppulses = "";
    if (!dumppulses.empty())
    {
        std::ofstream out(dumppulses.c_str());
        out << std::fixed;
        for (const Pulse& p : pulses)
            out << std::setprecision(6) << p.t + torg << " "
                << std::setprecision(3) << p.r(0) + rorg(0) << " "
                << p.r(1) + rorg(1) << " " << p.r(2) + rorg(2) << " " << p.d
                << " " << std::setprecision(6) << p.n(0) << " " << p.n(1) << " "
                << p.n(2) << "\n";
    }
    s = COMPLETE;
    Solve();
}

bool PulseCollection::EstimatedPositionVelocity(double t, Vector3d& r,
                                                Vector3d& v) const
{
    std::vector<Pulse> psub;
    for (const Pulse& p : pulses)
    {
        if (p.t >= t - m_tblock && p.t <= t + m_tblock)
        {
            psub.push_back(p);
            psub.back().t -= t;
        }
    }
    int k = int(psub.size());
    if (k < m_estn)
        return false;
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
    Matrix<double, Eigen::Dynamic, 6> A(2 * k, 6);
    VectorXd B(2 * k);
    for (int l = 0; l < k; ++l)
    {
        const Pulse& p = psub[l];
        A(2 * l + 0, 0) = p.n(2);
        A(2 * l + 0, 2) = -p.n(0);
        A(2 * l + 1, 1) = p.n(2);
        A(2 * l + 1, 2) = -p.n(1);
        A(2 * l + 0, 3) = p.n(2) * p.t;
        A(2 * l + 0, 5) = -p.n(0) * p.t;
        A(2 * l + 1, 4) = p.n(2) * p.t;
        A(2 * l + 1, 5) = -p.n(1) * p.t;
        A(2 * l + 0, 1) = A(2 * l + 0, 4) = A(2 * l + 1, 0) = A(2 * l + 1, 3) =
            0;
        B(2 * l + 0) = p.n(2) * p.r(0) - p.n(0) * p.r(2);
        B(2 * l + 1) = p.n(2) * p.r(1) - p.n(1) * p.r(2);
        A.row(2 * l + 0) *= p.d;
        B(2 * l + 0) *= p.d;
        A.row(2 * l + 1) *= p.d;
        B(2 * l + 1) *= p.d;
    }
    Matrix<double, 1, 6> rv(
        A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B));
    r = rv.head<3>();
    v = rv.tail<3>();
    return true;
}

Vector3d PulseCollection::EstimatedPosition(double t, double& estt) const
{
    int num = pulses.size();
    if (num == 0)
        throw std::runtime_error(
            "PulseCollection: no pulses for EstimatePosition");
    int k = -1;
    for (int i = 0; i < num; ++i)
    {
        if (pulses[i].t > t)
        {
            k = i;
            break;
        }
    }
    if (k < 0)
        k = num - 1;
    // equation for pulse starting at r in direction n, distance = s
    //
    //   x = rx + nx * s
    //   y = ry + ny * s
    //   z = rz + nz * s
    //
    // replace s by z as parameterization, s = (z - rz) / nz
    //
    //   nz * x         - nx * z = nz * rx - nx * rz
    //           nz * y - ny * z = nz * ry - ny * rz
    //  or
    //
    //    A . [x,y,z]' = B
    //
    //  where
    //
    //   A = [nz,  0, -nx]
    //       [ 0, nz, -ny]
    //       ... for all pulses
    //   B = [nz * rx - nx * rz]
    //       [nz * ry - ny * rz]
    //       ... for all pulses

    Matrix<double, Eigen::Dynamic, 3> A(2 * m_estn, 3);
    VectorXd B(2 * m_estn);
    int l = 0;
    double tsum = 0;
    for (int j = 0; j < num; ++j)
    {
        bool b = false;
        if (k - j >= 0)
        {
            b = true;
            const Pulse& p = pulses[k - j];
            A(2 * l + 0, 0) = p.n(2);
            A(2 * l + 0, 1) = 0;
            A(2 * l + 0, 2) = -p.n(0);
            A(2 * l + 1, 0) = 0;
            A(2 * l + 1, 1) = p.n(2);
            A(2 * l + 1, 2) = -p.n(1);
            B(2 * l + 0) = p.n(2) * p.r(0) - p.n(0) * p.r(2);
            B(2 * l + 1) = p.n(2) * p.r(1) - p.n(1) * p.r(2);
            tsum += p.t;
            ++l;
            if (l == m_estn)
                break;
        }
        if (j > 0 && k + j < num)
        {
            b = true;
            const Pulse& p = pulses[k + j];
            if (p.d > 0)
            {
                A(2 * l + 0, 0) = p.n(2);
                A(2 * l + 0, 1) = 0;
                A(2 * l + 0, 2) = -p.n(0);
                A(2 * l + 1, 0) = 0;
                A(2 * l + 1, 1) = p.n(2);
                A(2 * l + 1, 2) = -p.n(1);
                B(2 * l + 0) = p.n(2) * p.r(0) - p.n(0) * p.r(2);
                B(2 * l + 1) = p.n(2) * p.r(1) - p.n(1) * p.r(2);
                tsum += p.t;
                ++l;
                if (l == m_estn)
                    break;
            }
        }
        if (!b)
            break;
    }
    if (l < m_estn)
        throw std::runtime_error(
            "PulseCollection: insuffient multiple return pulses");

    estt = tsum / m_estn;
    Vector3d r(A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B));
    return r;
}

void PulseCollection::Solve()
{
    int npulses = int(pulses.size());
    if (npulses == 0)
        throw std::runtime_error("PulseCollection: no pulses for Solve");
    double tstart = std::floor(pulses[0].t / m_tblock);
    int num = int(std::ceil(pulses[npulses - 1].t / m_tblock) - tstart) - 1;
    if (num < 1)
        throw std::runtime_error("PulseCollection: no time interval for Solve");
    tstart *= m_tblock;
    traj = SplineFit3(num, m_tblock, tstart);
    std::string inittraj = "";
    if (inittraj.empty())
    {
        for (int i = 0; i <= num; ++i)
        {
            traj.missing[i] = !EstimatedPositionVelocity(tstart + i * m_tblock,
                                                         traj.r[i], traj.v[i]);
            traj.v[i] *= m_tblock;
        }
        if (!traj.fillmissing(true))
            throw std::runtime_error("PulseCollection: to few pulses for "
                                     "initial estimate of trajectory");
        std::string dumpinittraj = "";
        if (!dumpinittraj.empty())
        {
            std::ofstream str(dumpinittraj.c_str());
            str << std::fixed << std::setprecision(3);
            for (int i = 0; i <= num; ++i)
                str << tstart + i * m_tblock + torg << " "
                    << traj.r[i](0) + rorg(0) << " " << traj.r[i](1) + rorg(1)
                    << " " << traj.r[i](2) + rorg(2) << " "
                    << traj.v[i](0) / m_tblock << " " << traj.v[i](1) / m_tblock
                    << " " << traj.v[i](2) / m_tblock << " " << traj.missing[i]
                    << "\n";
        }
    }
    else
    {
        std::ifstream ftraj(inittraj.c_str());
        double t;
        Vector3d r, v;
        while (ftraj >> t >> r(0) >> r(1) >> r(2) >> v(0) >> v(1) >> v(2))
        {
            t -= torg;
            double f = (t - tstart) / m_tblock;
            int i = int(floor(f + 0.5));
            f -= i;
            if (abs(f) < 0.1 && i >= 0 && i <= num)
            {
                r -= rorg;
                v *= m_tblock;
                traj.r[i] = r;
                traj.v[i] = v;
            }
        }
    }

    google::InitGoogleLogging("LidarTrajectory");
    ceres::Problem problem;

    int skip = 1;
    int k = npulses / skip;
    // Set up residual blocks for pulses
    int n = 0;
    for (int j = 0; j < k; ++j)
    {
        const Pulse& p = pulses[j * skip];
        if (p.d == 0)
            continue;
        auto tconv = traj.tconvert(p.t);
        int i = tconv.first;
        double t = tconv.second;
        ceres::CostFunction* cost_function =
            new ceres::AutoDiffCostFunction<FirstLastError,
                                            2,          // number of residuals
                                            3, 3, 3, 3> // data for cubic fit
            (new FirstLastError(p, t));
        ceres::LossFunction* loss_function =
            (ceres::LossFunction*)(new ceres::CauchyLoss(m_dr));
        problem.AddResidualBlock(cost_function, loss_function, traj.r[i].data(),
                                 traj.v[i].data(), traj.r[i + 1].data(),
                                 traj.v[i + 1].data());
        ++n;
    }
    if (m_accelweight > 0)
    {
        for (int i = 1; i < num; ++i)
        {
            ceres::CostFunction* cost_function =
                new ceres::AutoDiffCostFunction<AccelJumpConstraint<3>,
                                                3, // number of residuals
                                                3, 3, 3, 3,
                                                3> // data for cubic fit
                (new AccelJumpConstraint<3>(m_tblock));
            ceres::LossFunction* loss_function =
                (ceres::LossFunction*)(new ceres::ScaledLoss(
                    NULL, m_accelweight, ceres::TAKE_OWNERSHIP));
            problem.AddResidualBlock(cost_function, loss_function,
                                     traj.r[i - 1].data(), traj.v[i - 1].data(),
                                     traj.v[i].data(), traj.r[i + 1].data(),
                                     traj.v[i + 1].data());
        }
    }
    if (m_clampweight > 0 || m_straddleweight > 0)
    {
        for (int i = 1; i < num; ++i)
        {
            double w = traj.missing[i] ? m_straddleweight : m_clampweight;
            if (!(w > 0))
                continue;
            ceres::CostFunction* cost_function =
                new ceres::AutoDiffCostFunction<ClampConstraint<3>,
                                                3, // number of residuals
                                                3, 3, 3, 3,
                                                3> // data for cubic fit
                (new ClampConstraint<3>(m_tblock));
            ceres::LossFunction* loss_function =
                (ceres::LossFunction*)(new ceres::ScaledLoss(
                    NULL, w, ceres::TAKE_OWNERSHIP));
            problem.AddResidualBlock(cost_function, loss_function,
                                     traj.r[i - 1].data(), traj.v[i - 1].data(),
                                     traj.r[i].data(), traj.r[i + 1].data(),
                                     traj.v[i + 1].data());
        }
    }
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    if (m_vlevel > 0)
        options.minimizer_progress_to_stdout = true;
    else
        options.logging_type = ceres::SILENT;
    options.max_linear_solver_iterations = m_niter;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    if (m_vlevel > 0)
        std::cerr << summary.FullReport() << "\n";
}

Vector3d PulseCollection::Trajectory(double t, Vector3d& v, Vector3d& a) const
{
    if (s != COMPLETE)
        throw std::runtime_error("PulseCollection: not yet finalized");
    Vector3d pos = traj.position(t - torg, v, a) + rorg;
    return pos;
}

} // namespace LidarTrajectory

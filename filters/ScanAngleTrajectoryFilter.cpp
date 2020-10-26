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

#include "ScanAngleTrajectoryFilter.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define deg2rad(d) ((d)*M_PI / 180.0)

namespace pdal
{

using namespace Dimension;

static StaticPluginInfo const s_info{
    "filters.scanangletrajectory", "Scan-angle trajectory estimation filter",
    "http://pdal.io/stages/filters.scanangletrajectory.html"};

CREATE_STATIC_STAGE(ScanAngleTrajectoryFilter, s_info)

ScanAngleTrajectoryFilter::ScanAngleTrajectoryFilter() {}

std::string ScanAngleTrajectoryFilter::getName() const
{
    return s_info.name;
}

void ScanAngleTrajectoryFilter::addArgs(ProgramArgs& args)
{
    // setPositional() Makes the argument required.
    args.add("trim_a", "Extent of extreme scan angles to remove (deg)",
             m_trim_a, 5);
    args.add("min_delta_a",
             "Minimum scan angle difference between point pairs (deg)",
             m_min_delta_a, 15);
    args.add("jitter",
             "Half-size of box filter applied to scan angles when identifying "
             "scan sweeps",
             m_jitter, 8);
    args.add("output_debug_files", "If true, output debug files.",
             m_output_debug_files, false);
    args.add("output_prefix", "Output filename prefix for debug files",
             m_output_prefix);
}

std::vector<int>
ScanAngleTrajectoryFilter::smoothScanAngles(const std::vector<int>& scan_angles,
                                            const int jitter)
{
    log()->get(LogLevel::Debug2)
        << "====== ScanAngleTrajectoryFilter::smoothScanAngles ====="
        << std::endl;

    int num_points = scan_angles.size();

    int filter_size = m_jitter * 2 + 1;

    log()->get(LogLevel::Debug2)
        << "Smoothing scan angle with filter size " << filter_size << std::endl;

    if (num_points < filter_size)
    {
        throwError("Error! number of points less than filter size.");
    }

    std::vector<int> a_smooth(num_points); // output of this function

    for (int i = jitter; i < num_points - jitter; i++)
    {
        // sum scan angles from data
        int sum = scan_angles[i];
        for (int j = 1; j <= jitter; j++)
        {
            sum += scan_angles[i + j];
            sum += scan_angles[i - j];
        }

        // average and round
        a_smooth[i] = round((float)sum / filter_size);
    }

    // fill in front and back
    int first_value = a_smooth[jitter];
    for (int i = 0; i < jitter; i++)
    {
        a_smooth[i] = first_value;
    }

    int last_value = a_smooth[num_points - jitter - 1];
    for (int i = num_points - jitter; i < num_points; i++)
    {
        a_smooth[i] = last_value;
    }

    log()->get(LogLevel::Debug2) << "Completed." << std::endl;

    return a_smooth;
}

std::vector<int> ScanAngleTrajectoryFilter::sweepIndices(
    const std::vector<int>& smooth_scan_angles)
{
    log()->get(LogLevel::Debug2)
        << "====== ScanAngleTrajectoryFilter::sweep_indices =====" << std::endl;

    int num_angles = smooth_scan_angles.size();
    log()->get(LogLevel::Debug4) << "num_angles = " << num_angles << std::endl;

    std::vector<int> d;
    std::vector<int> d2;
    std::vector<int> idx;
    std::vector<int> idx2;

    // d
    for (int i = 1; i < num_angles; i++)
    {
        int diff = smooth_scan_angles[i] - smooth_scan_angles[i - 1];
        if (diff > 0)
        {
            d.push_back(1);
            idx.push_back(i - 1);
        }
        else if (diff < 0)
        {
            d.push_back(-1);
            idx.push_back(i - 1);
        }
        else
        {
            // diff == 0. do nothing
            continue;
        }
    }

    log()->get(LogLevel::Debug4) << "num_da = " << d.size() << std::endl;

    // d2
    for (int i = 1; i < d.size(); i++)
    {
        int diff = d[i] - d[i - 1];
        if (diff > 0)
        {
            d2.push_back(1);
            idx2.push_back(i - 1);
        }
        else if (diff < 0)
        {
            d2.push_back(-1);
            idx2.push_back(i - 1);
        }
        else
        {
            // diff == 0. do nothing
            continue;
        }
    }

    log()->get(LogLevel::Debug4) << "num_d2a = " << d2.size() << std::endl;

    // indices
    std::vector<int> sweep_indices;
    for (int i = 0; i < idx2.size(); i++)
    {
        sweep_indices.push_back(idx[idx2[i] + 1] + 1);
    }

    log()->get(LogLevel::Debug2)
        << "Completed generating sweep_indices" << std::endl;
    log()->get(LogLevel::Debug2)
        << "Num sweep_indices = " << sweep_indices.size() << std::endl;

    if (m_output_debug_files)
    {
        // debug: write a in file
        std::ofstream outfile;
        outfile.open(m_output_prefix + "_sweep_indices.txt");
        for (int i = 0; i < sweep_indices.size(); i++)
        {
            outfile << sweep_indices[i] << std::endl;
            log()->get(LogLevel::Debug4) << sweep_indices[i] << std::endl;
        }

        outfile.close();
    }

    return sweep_indices;
}

bool ScanAngleTrajectoryFilter::getPointPairIndices(
    const std::vector<int> sweep_angles, int min_delta_a, int& low_idx1,
    int& low_idx2, int& high_idx1, int& high_idx2)
{
    int num_bins = 181;
    int min_angle = -90.5;
    int max_angle = 90.5;
    Histogram hist(sweep_angles, num_bins, min_angle, max_angle);
    std::vector<int> hist_vec = hist.get_hist();
    int hist_median = hist.get_median();

    std::vector<int> dense_bins;

    // dense bins
    for (int i = 0; i < hist_vec.size(); i++)
    {
        if (hist_vec[i] >= hist_median)
        {
            dense_bins.push_back(i);
        }
    }

    int min_dense_a = dense_bins[0] - 90;
    int max_dense_a = dense_bins[dense_bins.size() - 1] - 90;

    std::vector<int> low_idx;
    std::vector<int> high_idx;

    // check for sufficient scan angle geometry
    if (max_dense_a - min_dense_a >= min_delta_a)
    {
        for (int i = 0; i < sweep_angles.size(); i++)
        {
            // get indices where sweep_angles == min_dense_a
            if (sweep_angles[i] == min_dense_a)
            {
                low_idx.push_back(i);
            }
            // get indices where sweep_angles == max_dense_a
            if (sweep_angles[i] == max_dense_a)
            {
                high_idx.push_back(i);
            }
        }
    }

    if (low_idx.size() > 1 && high_idx.size() > 1)
    {
        low_idx1 = low_idx[0];
        low_idx2 = low_idx[low_idx.size() - 1];
        high_idx1 = high_idx[high_idx.size() - 1];
        high_idx2 = high_idx[0];

        return true;
    }

    return false;
}

ScanAngleTrajectoryFilter::Position
ScanAngleTrajectoryFilter::traj_xyz(const Position txyz_low, const int a_low,
                                    const Position txyz_high, const int a_high)
{
    double xl = txyz_low.x;
    double yl = txyz_low.y;
    double zl = txyz_low.z;
    int al = a_low;
    double xh = txyz_high.x;
    double yh = txyz_high.y;
    double zh = txyz_high.z;
    int ah = a_high;

    // L2-norm
    double b = std::sqrt((xh - xl) * (xh - xl) + (yh - yl) * (yh - yl) +
                         (zh - zl) * (zh - zl));
    double beta_l = deg2rad(90 + al);
    double beta_h = deg2rad(90 - ah);
    double alpha = deg2rad(ah - al);
    double gamma = beta_h - std::asin((zl - zh) / b);
    double d = (b * std::sin(gamma)) / std::sin(alpha);
    double theta = std::atan2(yh - yl, xh - xl);

    double x_est = xl + d * std::cos(beta_l) * std::cos(theta);
    double y_est = yl + d * std::cos(beta_l) * std::sin(theta);
    double z_est = zl + d * std::sin(beta_l);

    return Position(-1, x_est, y_est, z_est);
}

std::vector<ScanAngleTrajectoryFilter::Position>
ScanAngleTrajectoryFilter::computeTrajectory(
    const PointViewPtr view, const std::vector<int>& a,
    const std::vector<int>& sweep_indices, const int trim_a,
    const int min_delta_a)
{
    log()->get(LogLevel::Debug2)
        << "====== ScanAngleTrajectoryFilter::computeTrajectory ====="
        << std::endl;

    // debug files
    std::string filename = m_output_prefix + "_estimated_traj.txt";
    std::ofstream outfile;

    std::string filename2 = m_output_prefix + "_sweeps.txt";
    std::ofstream outfile2;

    if (m_output_debug_files)
    {
        outfile.open(filename);
        outfile2.open(filename2);
    }

    std::vector<Position> traj_txyz;

    int num_sweeps = sweep_indices.size();
    log()->get(LogLevel::Debug2) << "estimating trajectory from " << num_sweeps
                                 << " sweeps" << std::endl;
    log()->get(LogLevel::Debug2)
        << "trim_a = " << trim_a << ", min_delta_a = " << min_delta_a
        << std::endl;

    // for (int i = 0; i<num_sweeps - 1; i++)
    for (int i = 0; i < num_sweeps - 1; i++)
    {
        int idx_sweep_start = sweep_indices[i];
        int idx_sweep_end = sweep_indices[i + 1];

        // sweep data
        std::vector<int>::const_iterator start = a.begin() + idx_sweep_start;
        std::vector<int>::const_iterator end = a.begin() + idx_sweep_end;
        std::pair<std::vector<int>::const_iterator,
                  std::vector<int>::const_iterator>
            a_minmax = std::minmax_element(start, end);
        int a_min = *a_minmax.first;
        int a_max = *a_minmax.second;

        if (m_output_debug_files)
        {
            outfile2 << "sweep[" << i << "]: index = (" << idx_sweep_start
                     << ", " << idx_sweep_end - 1 << "), angle (min,max) = ("
                     << a_min << ", " << a_max << ")" << std::endl;
        }

        std::vector<int> a_sweep;
        std::vector<Position> txyz_sweep;

        // discard extreme scan angles
        for (int ia = idx_sweep_start; ia < idx_sweep_end; ia++)
        {
            if ((a_min + trim_a <= a[ia]) && (a[ia] <= a_max - trim_a))
            {
                a_sweep.push_back(a[ia]);

                double t = view->point(ia).getFieldAs<double>(Id::GpsTime);
                double x = view->point(ia).getFieldAs<double>(Id::X);
                double y = view->point(ia).getFieldAs<double>(Id::Y);
                double z = view->point(ia).getFieldAs<double>(Id::Z);

                txyz_sweep.push_back(Position(t, x, y, z));
            }
        }

        if (a_sweep.size() < 4)
        {
            // not sufficient data - skip to next sweep
            continue;
        }

        // get point pair indices
        int low_idx1 = -1;
        int low_idx2 = -1;
        int high_idx1 = -1;
        int high_idx2 = -1;
        if (getPointPairIndices(a_sweep, min_delta_a, low_idx1, low_idx2,
                                high_idx1, high_idx2))
        {
            log()->get(LogLevel::Debug3)
                << "point pair found (low, high): (" << low_idx1 << ", "
                << high_idx1 << "), (" << low_idx2 << ", " << high_idx2 << ")"
                << std::endl;

            if (m_output_debug_files)
            {
                outfile2 << "point pair found (low, high): (" << low_idx1
                         << ", " << high_idx1 << "), (" << low_idx2 << ", "
                         << high_idx2 << ")" << std::endl;
            }

            // mean trajectory solution in space and time
            Position p1_low = txyz_sweep[low_idx1];
            Position p1_high = txyz_sweep[high_idx1];

            int a1_low = a_sweep[low_idx1];
            int a1_high = a_sweep[high_idx1];

            Position txyz1_est = traj_xyz(p1_low, a1_low, p1_high, a1_high);

            // mean trajectory solution in space and time
            Position p2_low = txyz_sweep[low_idx2];
            Position p2_high = txyz_sweep[high_idx2];

            int a2_low = a_sweep[low_idx2];
            int a2_high = a_sweep[high_idx2];

            Position txyz2_est = traj_xyz(p2_low, a2_low, p2_high, a2_high);

            // mean
            double mean_x = (txyz1_est.x + txyz2_est.x) * 0.5;
            double mean_y = (txyz1_est.y + txyz2_est.y) * 0.5;
            double mean_z = (txyz1_est.z + txyz2_est.z) * 0.5;

            // time
            double mean_t =
                (txyz_sweep[low_idx1].t + txyz_sweep[low_idx2].t +
                 txyz_sweep[high_idx1].t + txyz_sweep[high_idx2].t) /
                4;

            traj_txyz.push_back(Position(mean_t, mean_x, mean_y, mean_z));

            if (m_output_debug_files)
            {
                outfile << std::fixed << std::setprecision(5) << mean_t << " "
                        << mean_x << " " << mean_y << " " << mean_z
                        << std::endl;
            }
        }
    }

    if (m_output_debug_files)
    {
        outfile.close();
        outfile2.close();
    }

    return traj_txyz;
}

PointViewSet ScanAngleTrajectoryFilter::run(const PointViewPtr view)
{
    log()->get(LogLevel::Debug2)
        << "====== ScanAngleTrajectoryFilter::run =====" << std::endl;
    log()->get(LogLevel::Debug2) << "traim_a: " << m_trim_a << std::endl;
    log()->get(LogLevel::Debug2)
        << "min_delta_a: " << m_min_delta_a << std::endl;
    log()->get(LogLevel::Debug2) << "jitter :" << m_jitter << std::endl;

    if (m_output_debug_files && m_output_prefix == "")
    {
        m_output_prefix = "debug_output";
    }

    // scan angle rank
    // int num_points = view->size();
    int num_points = 20000;
    m_scan_angles = std::vector<int>(num_points);
    for (int i = 0; i < num_points; i++)
    {
        m_scan_angles[i] = view->point(i).getFieldAs<int>(Id::ScanAngleRank);
    }

    m_smooth_scan_angles = smoothScanAngles(m_scan_angles, m_jitter);
    m_sweep_indices = sweepIndices(m_smooth_scan_angles);
    m_txyz_source_est = computeTrajectory(view, m_scan_angles, m_sweep_indices,
                                          m_trim_a, m_min_delta_a);

    log()->get(LogLevel::Debug2)
        << "trajectory estimation completed" << std::endl;

    // outViewSet
    m_outViewSet.clear();
    PointViewPtr new_view = view->makeNew();
    for (std::vector<Position>::iterator ptr = m_txyz_source_est.begin();
         ptr < m_txyz_source_est.end(); ptr++)
    {
        PointId cnt = new_view->size();
        new_view->setField(Id::GpsTime, cnt, ptr->t);
        new_view->setField(Id::X, cnt, ptr->x);
        new_view->setField(Id::Y, cnt, ptr->y);
        new_view->setField(Id::Z, cnt, ptr->z);
    }
    m_outViewSet.insert(new_view);
    log()->get(LogLevel::Debug2)
        << "created new view with size " << new_view->size() << std::endl;

    return m_outViewSet;
}

} // namespace pdal

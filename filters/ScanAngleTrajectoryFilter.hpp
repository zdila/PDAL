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

#include <pdal/Filter.hpp>

namespace pdal
{
class PDAL_DLL ScanAngleTrajectoryFilter : public Filter
{
    struct Position
    {
        double t;
        double x;
        double y;
        double z;

        Position(double t_, double x_, double y_, double z_)
        {
            t = t_;
            x = x_;
            y = y_;
            z = z_;
        };
    };

    class Histogram
    {
        std::vector<int> m_hist;
        double m_dt;
        double m_min;
        double m_max;

    public:
        Histogram(const int num_buckets, const double min_value,
                  const double max_value)
        {
            m_hist = std::vector<int>(num_buckets, 0);
            m_min = min_value;
            m_max = max_value;
            m_dt = (m_max - m_min) / num_buckets;
        }

        template <class T>
        Histogram(const std::vector<T> vec, const int num_buckets,
                  const double min_value, const double max_value)
        {
            m_hist = std::vector<int>(num_buckets, 0);
            m_min = min_value;
            m_max = max_value;
            m_dt = (m_max - m_min) / num_buckets;

            int vec_length = vec.size();
            for (int i = 0; i < vec_length; i++)
            {
                insert_value(vec[i]);
            }
        }

        void insert_value(const double val)
        {
            int index = (val - m_min) / m_dt;
            if (0 <= index && index < m_hist.size())
            {
                m_hist[index] += 1;
            }
        }

        std::vector<int> get_hist()
        {
            return m_hist;
        }

        int get_median()
        {
            std::vector<int> non_zero_hist;
            // num-zero items
            for (int i = 0; i < m_hist.size(); i++)
            {
                if (m_hist[i] > 0)
                {
                    non_zero_hist.push_back(m_hist[i]);
                }
            }

            std::vector<int>::iterator b = non_zero_hist.begin();
            std::vector<int>::iterator e = non_zero_hist.end();
            std::vector<int>::iterator med = b;
            std::advance(med, non_zero_hist.size() / 2);

            // This makes the 2nd position hold the median.
            std::nth_element(b, med, e);

            int median = *med;
            return median;
        }
    };

public:
    ScanAngleTrajectoryFilter();
    ScanAngleTrajectoryFilter&
    operator=(const ScanAngleTrajectoryFilter&) = delete;
    ScanAngleTrajectoryFilter(const ScanAngleTrajectoryFilter&) = delete;

    std::string getName() const;

private:
    SpatialReference m_srs;
    PointViewSet m_outViewSet;

    // args
    std::string m_output_prefix;
    int m_trim_a;      // Extent of extreme scan angles to remove (deg)
    int m_min_delta_a; // Minimum scan angle difference between point pairs
                       // (deg)
    int m_jitter;      // Half-size of box filter applied to scan angles when
                       // identifying scan sweeps
    bool m_output_debug_files;

    std::vector<int> m_scan_angles;
    std::vector<int> m_smooth_scan_angles;
    std::vector<int> m_sweep_indices;

    // estimated result
    std::vector<Position> m_txyz_source_est;

    virtual void addArgs(ProgramArgs& args);
    virtual PointViewSet run(const PointViewPtr view);

    std::vector<int> smoothScanAngles(const std::vector<int>& scan_angles,
                                      const int jitter);
    std::vector<int> sweepIndices(const std::vector<int>& smooth_scan_angles);
    std::vector<Position>
    computeTrajectory(const PointViewPtr view, const std::vector<int>& a,
                      const std::vector<int>& sweep_indices, const int trim_a,
                      int const min_delta_a);
    bool getPointPairIndices(const std::vector<int> sweep_angles, int minDeltaA,
                             int& low_idx1, int& low_idx2, int& high_idx1,
                             int& high_idx2);
    Position traj_xyz(const Position txyz_low, const int a_low,
                      const Position txyz_high, const int a_high);
};

} // namespace pdal

/******************************************************************************
 * Copyright (c) 2021, Bradley J Chambers (brad.chambers@gmail.com)
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

// PDAL implementation of subsampling strategy presented in Groh F.,
// Wieschollek P., Lensch H.P.A. (2019) Flex-Convolution. In: Jawahar C., Li
// H., Mori G., Schindler K. (eds) Computer Vision – ACCV 2018. ACCV 2018.
// Lecture Notes in Computer Science, vol 11361. Springer, Cham.
// https://doi.org/10.1007/978-3-030-20887-5_7

#include "InverseDensityImportanceSamplingFilter.hpp"

#include "private/Segmentation.hpp"

#include <pdal/KDIndex.hpp>
#include <pdal/util/ProgramArgs.hpp>

#include <algorithm>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace pdal
{

static PluginInfo const s_info{"filters.idis",
                               "Inverse density importance sampling filter",
                               "http://pdal.io/stages/filters.idis.html"};

CREATE_STATIC_STAGE(InverseDensityImportanceSamplingFilter, s_info)

std::string InverseDensityImportanceSamplingFilter::getName() const
{
    return s_info.name;
}

InverseDensityImportanceSamplingFilter::InverseDensityImportanceSamplingFilter()
{
}

void InverseDensityImportanceSamplingFilter::addArgs(ProgramArgs& args)
{
    args.add("count", "Target number of points after sampling", m_count,
             point_count_t(1000));
    args.add("knn", "Number of neighbors to consider while computing density",
             m_knn, point_count_t(16));
}

PointViewSet InverseDensityImportanceSamplingFilter::run(PointViewPtr inView)
{
    // Return empty PointViewSet if the input PointView has no points.
    PointViewSet viewSet;
    if (!inView->size() || (inView->size() < m_count))
        return viewSet;

    PointIdList ids =
        Segmentation::inverseDensityImportanceSampling(*inView, m_count, m_knn);

    PointViewPtr outView = inView->makeNew();
    for (auto const& id : ids)
        outView->appendPoint(*inView, id);
    viewSet.insert(outView);
    return viewSet;
}

} // namespace pdal

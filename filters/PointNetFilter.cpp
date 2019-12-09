/******************************************************************************
 * Copyright (c) 2020, Bradley J Chambers (brad.chambers@gmail.com)
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

#include "PointNetFilter.hpp"

#include <pdal/EigenUtils.hpp>
#include <pdal/KDIndex.hpp>
#include <pdal/util/FileUtils.hpp>
#include <pdal/util/ProgramArgs.hpp>

#include "private/DimRange.hpp"
#include "private/Segmentation.hpp"

#include <Eigen/Dense>

#include <torch/torch.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace pdal
{

using namespace Dimension;
using namespace Eigen;

static StaticPluginInfo const s_info{
    "filters.pointnet2", "PointNet++ (Qi et al., 201x)",
    "http://pdal.io/stages/filters.pointnet2.html"};

CREATE_STATIC_STAGE(PointNetFilter, s_info)

struct PointNetArgs
{
};

struct SAModule : torch::nn::Module
{
	// ratio
	// r
	// conv = PointConv(nn)
	
	torch::Tensor forward(torch::Tensor x)
	{
		// idx = fps(pos, batch, ratio=ratio)
		// row, col = radius(pos, pos[idx], r, batch, batch[idx], max_num_neighbors=64)
		// edge_index = torch.stack([col, row], dim=0)
		// x = conv(x, (pos, pos[idx]), edge_index)
		// pos, batch = pos[idx], batch[idx]
		// return x, pos, batch
	}
};

struct GlobalSAModule : torch::nn::Module
{
	// nn
	
	torch::Tensor forward(torch::Tensor x)
	{
		// x = nn(torch.cat([x, pos], dim=1))
		// x = global_max_pool(x, batch)
		// pos = pos.new_zeros((x.size(0), 3))
		// batch = torch.arange(x.size(0), device=batch.device)
		// return x, pos, batch
	}
};

torch::nn::Sequential MLP(
		// torch::nn::Linear(channels[i-1], channels[i]),
		// torch::nn::Functional(torch::relu),
		// torch::BatchNorm(i)
		);
//MLP->to(device);

struct Net : torch::nn::Module
{
    // sa1_module = SAModule(0.2, 0.2, MLP([3, 64, 64, 128]))
    // sa2_module = SAModule(0.25, 0.4, MLP([128+3, 128, 128, 256]))
    // sa3_module = GlobalSAModule(MLP([256+3, 256, 512, 1024]))
    // fp3_module = FPModule(1, MLP([1024+256, 256, 256]))
    // fp2_module = FPModule(3, MLP([256+128, 256, 128]))
    // fp1_module = FPModule(3, MLP([128, 128, 128, 128]))
    // lin1 = torch.nn.Linear(128, 128)
    // lin2 = torch.nn.Linear(128, 128)
    // lin3 = torch.nn.Linear(128, num_classes)
    Net() : lin1(128, 128), lin2(128, 128), lin3(128, 6)
    {
	    // register_model("lin1", lin1);
	    // register_model("lin2", lin2);
	    // register_model("lin3", lin3);
    }

    // Forward
    torch::Tensor forward(torch::Tensor x)
    {
	    // sa0_out = (data.x, data.pos, data.batch)
	    // sa1_out = sa1_module(*sa0_out)
	    // sa2_out = sa2_module(*sa1_out)
	    // sa3_out = sa3_module(*sa2_out)
	    // fp3_out = fp3_module(*sa3_out, *sa2_out)
	    // fp2_out = fp2_module(*fp3_out, *sa1_out)
	    // x = fp1_module(*fp2_out, *sa0_out)
	    x = torch::relu(lin1(x));
	    x = torch::dropout(x, /*p=*/0.5, /*training=*/is_training());
	    x = lin2(x);
	    x = torch::dropout(x, /*p=*/0.5, /*training=*/is_training());
	    x = lin3(x);
	    return torch::log_softmax(x, /*dim=*/1); // maybe it should be -1, read docs
    }

    torch::nn::Linear lin1;
    torch::nn::Linear lin2;
    torch::nn::Linear lin3;
};

PointNetFilter::PointNetFilter() : m_args(new PointNetArgs) {}

PointNetFilter::~PointNetFilter() {}

std::string PointNetFilter::getName() const
{
    return s_info.name;
}

void PointNetFilter::addArgs(ProgramArgs& args)
{
}

void PointNetFilter::addDimensions(PointLayoutPtr layout)
{
    layout->registerDim(Id::Classification);
}

void PointNetFilter::prepared(PointTableRef table)
{
    const PointLayoutPtr layout(table.layout());
}

void PointNetFilter::ready(PointTableRef table)
{
}

PointViewSet PointNetFilter::run(PointViewPtr view)
{
    PointViewSet viewSet;
    if (!view->size())
        return viewSet;

    return viewSet;
}

} // namespace pdal

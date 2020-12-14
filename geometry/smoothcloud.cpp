/**
 * Copyright (c) 2017 Melown Technologies SE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * *  Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * *  Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/**
 *  @file smoothcloud.cpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Cloud smoothing stuff.
 *
 *  2014-07-09 (vasek) Stolen from window-mesh
 */

#include "utility/openmp.hpp"
#include "utility/progress.hpp"

#include "math/geometry_core.hpp"

#include "neighbors.hpp"
#include "smoothcloud.hpp"

namespace geometry {

math::Points3 CloudSmoother::operator()(const math::Points3 &points
                                        , const KdTree<math::Point3> &kdtree)
    const
{
    // sanity check
    if (!params_.feasible(points.size())) {
        return points;
    }

    LOG(info3) << "Smoothing cloud.";
    KdTree<math::Point3>::Neighbors neighbors;

    auto size(points.size());

    // allocate smooth cloud
    math::Points3 smooth(size);

    utility::ts::Progress progress("Cloud smoothing", size);

    UTILITY_OMP(parallel for private(neighbors))
    for (int64_t i = 0; i < size; ++i) {
        const auto &point(points[i]);
        auto &spoint(smooth[i]);

        // find neighbors
        collectNeighbors(kdtree, point, neighbors, params_.neighbors
                         , params_.radius, true);

        // calculate smoothed point
        math::Point3 sum;
        double wsum = 0.0;
        for (const auto &pt : neighbors) {
            sum += pt.first;
            wsum += 1.0;
        }
        spoint = sum * (1.0 / wsum);

        ++progress;
    }

    return smooth;
}

} // namespace geometry

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

#include "./neighbors.hpp"
#include "./smoothcloud.hpp"

namespace geometry {

math::Points3 CloudSmoother::operator()(const math::Points3 &points
                                        , const KdTree<math::Point3> &kdtree)
    const
{
    LOG(info3) << "Smoothing cloud.";
    KdTree<math::Point3>::Neighbors neighbors;

    auto size(points.size());

    if (size < params_.neighbors) {
        LOG(warn2) << "Too few points to smooth cloud. Keeping intact.";
        return points;
    }

    // allocate smooth cloud
    math::Points3 smooth(size);

    utility::ts::Progress progress("Cloud smoothing", size);

    UTILITY_OMP(parallel for private(neighbors))
    for (std::size_t i = 0; i < size; ++i) {
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

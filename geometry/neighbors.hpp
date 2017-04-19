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
 *  @file neighbors.hpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Neighbor search.
 *
 *  2014-07-08 (vasek) Stolen from window-mesh
 */

#ifndef geometry_neighbors_hpp_included_
#define geometry_neighbors_hpp_included_

#include "dbglog/dbglog.hpp"
#include "math/geometry_core.hpp"

#include "./kdtree.hpp"

namespace geometry {

/** Find neighbors of a point in given radius.
 */
template <typename PointType>
double collectNeighbors(const KdTree<PointType> &kdtree
                        , const PointType &point
                        , typename KdTree<PointType>::Neighbors &neighbors
                        , size_t max, double radius
                        , bool dontCutFirstRadius);

// implementation

template <typename PointType>
inline double
collectNeighbors(const KdTree<PointType> &kdtree
                 , const PointType &point
                 , typename KdTree<PointType>::Neighbors &neighbors
                 , size_t max, double radius, bool dontCutFirstRadius)
{
    typedef typename KdTree<PointType>::Neighbor Neighbor;

    // we need to limit max number of neighbors to same sane value
    max = std::min(max, kdtree.size() / 2);

    int iterations = 0;
    do {
        LOG(info1) << "collectNeighbors: using radius = " << radius;
        neighbors.clear();
        kdtree.template range<false>(point, radius, neighbors);
        iterations++;
        radius *= 2;
        LOG(info1) << "collectNeighbors: found = " << neighbors.size();
    } while (neighbors.size() < max);

    LOG(info1) << "! collectNeighbors: found = " << neighbors.size();

    if ((neighbors.size() > max) && (!dontCutFirstRadius || iterations > 1))
    {
        LOG(info1) << "too many points -> sorting and cutting";
        // sort and cut
        std::nth_element(neighbors.begin(), neighbors.begin() + max
                         , neighbors.end()
                         , [](const Neighbor &l, const Neighbor &r)
                         {
                             return l.second < r.second;
                         });
        neighbors.resize(max);
    }

    // make radius 20 % bigger for next round
    radius = std::sqrt(neighbors.rbegin()->second) * 1.2;
    LOG(info1) << "new radius = " << radius;
    return radius;
}

} // namespace geometry

#endif //geometry_neighbors_hpp_included_

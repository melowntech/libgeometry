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

/** Neighbor of a point: neighbor and its distance from point
 */
template <typename PointType>
struct Neighbor : std::pair<PointType, double>;

/** List of neighbors.
 */

template <typename PointType>
struct Neighbors : template std::vector<Neighbor<PointType> >;

/** Find neighbors of a point in given radius.
 */
template <typename PointType>
double collectNeighbors(const KdTree<PointType> &kdtree
                        , const PointType &point
                        , Neighbors<PointType> &neighbors, const size_t max
                        , double radius, bool dontCutFirstRadius);

// implementation

inline double collectNeighbors(const KdTree<PointType> &kdtree
                               , const PointType &point
                               , Neighbors<PointType> &neighbors
                               , const size_t max
                               , double radius, bool dontCutFirstRadius)
{
    int iterations = 0;
    do
    {
        LOG(info1) << "collectNeighbors: using radius = " << radius;
        neighbors.clear();
        kdtree.range<false>(point, radius, neighbors);
        iterations++;
        radius *= 2;
        LOG(info1) << "collectNeighbors: found = " << neighbors.size();
    }
    while (neighbors.size() < max);

    LOG(info1) << "! collectNeighbors: found = " << neighbors.size();

    if ((neighbors.size() > max) && (!dontCutFirstRadius || iterations > 1))
    {
        LOG(info1) << "too many points -> sorting and cutting";
        // sort and cut
        std::nth_element(neighbors.begin(), neighbors.begin() + max
                         , neighbors.end()
                         , [](const Neighbor<PointType> &l
                              , const Neighbor<PointType> &r)
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

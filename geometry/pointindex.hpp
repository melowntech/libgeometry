/**
 * @file pointindex.hpp
 * @author Jakub Cerveny <jakub.cerveny@ext.citationtech.net>
 *
 * Helper class for removing duplicate vertices.
 */

#ifndef geometry_pointindex_hpp_included_
#define geometry_pointindex_hpp_included_

#include <map>

namespace geometry {

//! Helper class for removing duplicate vertices.
//!
template<typename Point>
class PointIndex
{
    std::map<Point, unsigned> map_;
    unsigned next_;

public:
    PointIndex() : next_(0) {}

    //! Assigns (and returns) a unique index for the given point.
    //! Idices are reused for points that have been seen before.
    unsigned assign(const Point& pt)
    {
        auto pair = map_.insert(std::make_pair(pt, next_));
        if (pair.second) next_++;
        return pair.first->second;
    }
};

} // namespace geometry

#endif // geometry_pointindex_hpp_included_

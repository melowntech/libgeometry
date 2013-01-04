/**
 *  @file polyclip.hpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Polygon clipping stuff
 *
 *  2013-01-03 (vasek)    created
 */

#ifndef geometry_kdtree_hpp_included_
#define geometry_kdtree_hpp_included_

#include "math/geometry_core.hpp"

namespace geometry {

namespace detail {
    math::Points2 clip(const math::Viewport2f &viewport
                       , const math::Points2 &polygon);
}

template <typename T, typename U>
std::vector<math::Point2_<T> >
clip(const math::Viewport2_<U> &viewport
     , const std::vector<math::Point2_<T> > &polygon)
{
    auto res(detail::clip(math::Viewport2f(viewport)
                          , math::Points2(polygon.begin(), polygon.end())));
    return { res.begin(), res.end() };
}

} // namespace geometry

#endif // geometry_kdtree_hpp_included_

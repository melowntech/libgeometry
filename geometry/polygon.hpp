/**
 *  @file polyclip.hpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Polygon clipping stuff
 *
 *  2013-01-03 (vasek)    created
 */

#ifndef geometry_polygon_hpp_included_
#define geometry_polygon_hpp_included_

#include "math/geometry.hpp"

namespace geometry {

/** Clip polygon by provided viewport.
 *
 *  Returns polygon of same type as is input. Internally works with double
 *  precision, though.
 */
template <typename T, typename U>
std::vector<math::Point2_<T> >
clip(const math::Viewport2_<U> &viewport
     , const std::vector<math::Point2_<T> > &polygon);


/** Determine whether polygon is convex.
 */
template <typename T>
bool convex(const std::vector<math::Point2_<T> > &polygon);

// implementation

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

namespace detail {

inline char getSign(double value) {
    if (value > 0.) {
        return 1;
    } else if (value < 0.) {
        return 2;
    } else {
        return 0;
    }
}

template <typename T>
inline double cp(const math::Point2_<T> &prev
                 , const math::Point2_<T> &cur
                 , const math::Point2_<T> &next)
{
    return math::crossProduct(math::Point2_<T>(cur - prev)
                              , math::Point2_<T>(next - cur));
}

} // namespace detail

template <typename T>
bool convex(const std::vector<math::Point2_<T> > &polygon)
{
    if (polygon.size() < 3) {
        // this is not a polygon
        return false;
    }

    char sign(0);

    // first, special cases
    sign = detail::getSign(detail::cp(polygon.back(), polygon.front()
                                       , polygon[1]));

    sign |= detail::getSign(detail::cp(polygon[polygon.size() - 2]
                                       , polygon.back()
                                       , polygon.front()));
    // check for signs
    if (sign == 3) {
        return false;
    }

    // rest of polygon
    for (auto ip(polygon.begin() + 1), ep(polygon.end() - 1); ip != ep; ++ip) {
        sign |= detail::getSign(detail::cp(*(ip - 1), *ip, *(ip + 1)));
        if (sign == 3) {
            return false;
        }
    }

    return true;
}

} // namespace geometry

#endif // geometry_polygon_hpp_included_

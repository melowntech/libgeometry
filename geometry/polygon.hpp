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

/** Compute polygon area (via shoelace formula)
 */
template <typename T>
double area( const std::vector<math::Point2_<T> > & polygon  );

/** Clip polygon by provided viewport.
 *
 *  Returns polygon of same type as is input. Internally works with double
 *  precision, though.
 */
template <typename T, typename U>
std::vector<math::Point2_<T> >
clip(const math::Viewport2_<U> &viewport
     , const std::vector<math::Point2_<T> > &polygon);

/** Clip polygon by provided viewport.
 *
 *  Returns vector of same points as has input polygon. Internally works with
 *  double precision, though.
 */
template <typename Iterator, typename U>
std::vector<typename std::iterator_traits<Iterator>::value_type>
clip(const math::Viewport2_<U> &viewport, Iterator begin, Iterator end);


/** Determine whether polygon is convex.
 */
template <typename T>
bool convex(const std::vector<math::Point2_<T> > &polygon);

/** Determine whether point is inside CONVEX polygon.
 */
template <typename T>
bool inside(const std::vector<math::Point2_<T> > &polygon
            , const math::Point2_<T> &point);

/** Determine whether point is inside CONVEX polygon.
 */
template <typename PointType1, typename PointType2>
bool insidePolygon(const std::vector<PointType1> &polygon
                   , const PointType2 &point);

/****** implementation ******/

template <typename T>
double area( const std::vector<math::Point2_<T> > & polygon  ) {

    auto n = polygon.size();
    double retval( 0.0 );

    for ( decltype(n) i = 0; i < n; i++ ) {

        retval += polygon[i](0) * (
            polygon[ ( i + 1 ) % n ](1)
            - polygon[ ( i + n - 1 ) % n ](1) );
    }

    return 0.5 * retval;
}


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

template <typename Iterator, typename U>
std::vector<typename std::iterator_traits<Iterator>::value_type>
clip(const math::Viewport2_<U> &viewport
     , Iterator begin, Iterator end)
{
    auto res(detail::clip(math::Viewport2f(viewport)
                          , math::Points2(begin, end)));
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

template <typename T>
inline double cp(const math::Point3_<T> &prev
                 , const math::Point3_<T> &cur
                 , const math::Point3_<T> &next)
{
    return math::crossProduct(math::Point2_<T>(cur - prev)
                              , math::Point2_<T>(next - cur));
}

} // namespace detail

template <typename T>
inline bool convex(const std::vector<math::Point2_<T> > &polygon)
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

template <typename T>
inline bool inside(const std::vector<math::Point2_<T> > &polygon
                   , const math::Point2_<T> &point)
{
    // NB: works only for convex polygons
    if (polygon.size() < 3) {
        // this is not a polygon
        return false;
    }

    char sign(0);

    // first, special cases
    sign = detail::getSign(detail::cp(polygon.back(), polygon.front()
                                       , point));

    // rest of polygon
    for (auto ip(polygon.begin() + 1), ep(polygon.end()); ip != ep; ++ip) {
        sign |= detail::getSign(detail::cp(*(ip - 1), *ip, point));
        if (sign == 3) {
            return false;
        }
    }

    return true;
}

template <typename PolygonType, typename PointType2>
bool insidePolygon(const PolygonType &polygon
                   , const PointType2 &point)
{
    // NB: works only for convex polygons
    if (polygon.size() < 3) {
        // this is not a polygon
        return false;
    }

    char sign(0);

    // first, special cases
    auto ip(polygon.begin());
    sign = detail::getSign(detail::cp(polygon.back(), *ip, point));

    // rest of polygon
    auto pip(ip);
    ++ip;

    for (auto ep(polygon.end()); ip != ep; ++ip, ++pip) {
        sign |= detail::getSign(detail::cp(*pip, *ip, point));
        if (sign == 3) {
            return false;
        }
    }

    return true;
}

} // namespace geometry

#endif // geometry_polygon_hpp_included_

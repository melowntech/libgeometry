/**
 * Copyright (c) 2023 Melown Technologies SE
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
 *  @file polyop.cpp
 *  @author Jakub Cerveny <jakub.cerveny@melowntech.com>
 *
 *  Polygon operations, based on Boost.Geometry.
 *
 */

#include "polyop.hpp"
#include "polygon.hpp"
#include "boost-geometry-convert.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

namespace bg = boost::geometry;

namespace geometry {

typedef bg::model::d2::point_xy<double> bgPoint;
typedef bg::model::polygon<bgPoint, false, false> bgPolygon; // ccw, unclosed
typedef bg::model::multi_polygon<bgPolygon> bgMultiPolygon;
typedef bgPolygon::ring_type bgRing;

math::MultiPolygon intersectPolygons(const math::MultiPolygon &mp1,
                                     const math::MultiPolygon &mp2)
{
    bgMultiPolygon out;
    bg::intersection(convert2bg<bgMultiPolygon>(mp1), convert2bg<bgMultiPolygon>(mp2), out);
    return convert2math(out);
}

math::MultiPolygon subtractPolygons(const math::MultiPolygon &mp1,
                                    const math::MultiPolygon &mp2)
{
    bgMultiPolygon out;
    bg::difference(convert2bg<bgMultiPolygon>(mp1), convert2bg<bgMultiPolygon>(mp2), out);
    return convert2math(out);
}

math::MultiPolygon unitePolygons(const math::MultiPolygon &mp1,
                                 const math::MultiPolygon &mp2)
{
    bgMultiPolygon out;
    bg::union_(convert2bg<bgMultiPolygon>(mp1), convert2bg<bgMultiPolygon>(mp2), out);
    return convert2math(out);
}

math::MultiPolygon offsetPolygon(const math::MultiPolygon &mpoly,
                                 double distance, double miterLimit,
                                 double epsSimplify)
{
    bg::strategy::buffer::distance_symmetric<double> s_distance(distance);
    bg::strategy::buffer::side_straight s_side;
    bg::strategy::buffer::join_miter s_join(miterLimit);
    bg::strategy::buffer::end_flat s_end;
    bg::strategy::buffer::point_circle s_circle(1);

    bgMultiPolygon out;
    bg::buffer(convert2bg<bgMultiPolygon>(mpoly), out,
               s_distance, s_side, s_join, s_end, s_circle);

    if (epsSimplify > 0.0)
    {
        // get rid of collinear points produced by Boost.Geomertry around joins
        bgMultiPolygon out2;
        bg::simplify(out, out2, epsSimplify);
        return convert2math(out2);
    }
    else
    {
        return convert2math(out);
    }
}

} // namespace geometry

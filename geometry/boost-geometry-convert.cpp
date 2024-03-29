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
 *  @file boost-geometry-convert.cpp
 *
 *  Conversions from and to Boost.Geometry.
 *
 */

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

template<>
bgMultiPolygon convert2bg<bgMultiPolygon>(const math::MultiPolygon &mpoly)
{
    auto makeRing = [](const math::Points2d &points)
    {
        bgRing ring;
        ring.reserve(points.size());
        for (const auto &p : points) {
            ring.emplace_back(p(0), p(1));
        }
        return ring;
    };

    bgMultiPolygon result;
    result.reserve(mpoly.size());

    // add polygons with positive area as regular bgPolygons
    for (const auto &poly : mpoly) {
        if (geometry::area(poly) > 0) {
            result.push_back(bgPolygon({makeRing(poly)}));
        }
    }

    // add all holes as interior rings of the first bgPolygon
    // -- boost::geometry seems OK with that
    for (const auto &poly : mpoly) {
        if (geometry::area(poly) < 0) {
            bg::interior_rings(result.front()).push_back(makeRing(poly));
        }
    }

    return result;
}

template<>
math::MultiPolygon convert2math<bgMultiPolygon>(const bgMultiPolygon &bgmpoly)
{
    auto ringPoints = [](const bgPolygon::ring_type &ring)
    {
        math::Points2d result;
        result.reserve(ring.size());
        for (const auto &p : ring) {
            result.emplace_back(p.x(), p.y());
        }
        return result;
    };

    math::MultiPolygon result;
    result.reserve(bgmpoly.size());
    for (const auto &poly : bgmpoly)
    {
        result.push_back(ringPoints(poly.outer()));
        for (const auto &ring : poly.inners()) {
            result.push_back(ringPoints(ring));
        }
    }
    return result;
}

} // namespace geometry

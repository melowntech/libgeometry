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

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#if BOOST_VERSION < 105600
#  include <boost/geometry/multi/geometries/multi_polygon.hpp>
#else
#  include <boost/geometry/geometries/multi_polygon.hpp>
#endif

#include "dbglog/dbglog.hpp"

#include "nonconvexclip.hpp"
#include "triangulate.hpp"
#include "boost-geometry-convert.hpp"

namespace bg = boost::geometry;

namespace geometry {

typedef bg::model::d2::point_xy<double> Point;
typedef bg::model::polygon<Point, false, false> Polygon; // ccw, unclosed
typedef bg::model::multi_polygon<Polygon> MultiPolygon;
typedef Polygon::ring_type Ring;

namespace {

inline math::Points2d ringPoints(const Ring &ring)
{
    math::Points2d result;
    result.reserve(ring.size());
    for (const auto &p : ring) {
        result.emplace_back(p.x(), p.y());
    }
    return result;
}

inline double checkCcw(const math::Point2 &a, const math::Point2 &b
                       , const math::Point2 &c)
{
    return math::crossProduct(math::Point2(math::normalize(b - a))
                              , math::Point2(math::normalize(c - a)));
}

inline double area(const math::Triangle2d &t) {
    return 0.5 * std::abs(math::crossProduct( math::Point2(t[1] - t[0])
                                            , math::Point2(t[2] - t[0])));
}

} // namespace

math::Triangles3d clipTriangleNonconvex(const math::Triangle3d &tri_,
                                        const math::MultiPolygon &clipRegion)
{
    math::Triangle3d tri(tri_);

    // tri -> tri2
    math::Triangle2d tri2;
    for (int i = 0; i < 3; i++) {
        tri2[i](0) = tri[i](0);
        tri2[i](1) = tri[i](1);
    }

    bool flip = false;
    double ccw(checkCcw(tri2[0], tri2[1], tri2[2]));

    // convert clipRegion to boost MultiPolygon
    MultiPolygon clipMultiPoly { convert2bg(clipRegion) };

    if (std::abs(ccw) < 1e-4)
    {
        // TODO: handle exactly vertical triangles properly
        if (bg::within(Point{tri[0](0), tri[0](1)}, clipMultiPoly) &&
            bg::within(Point{tri[1](0), tri[1](1)}, clipMultiPoly) &&
            bg::within(Point{tri[2](0), tri[2](1)}, clipMultiPoly))
        {
            LOG(debug)
                << "Including near vertical triangle, all vertices lie inside.";
            return {tri};
        } else {
            LOG(debug) << "Excluding near vertical triangle, "
                          "not all vertices lie inside.";
            return {};
        }
    }

    // ensure counter-clockwise orientation for clipping
    if (ccw < 0.0) {
        std::swap(tri[1], tri[2]);
        std::swap(tri2[1], tri2[2]);
        flip = true;
    }

    // convert input to 2D polygons
    Polygon trianglePoly;
    bg::append(trianglePoly.outer(), Point(tri[0][0], tri[0][1]));
    bg::append(trianglePoly.outer(), Point(tri[1][0], tri[1][1]));
    bg::append(trianglePoly.outer(), Point(tri[2][0], tri[2][1]));

    // calculate intersection
    std::deque<Polygon> isect;
    bg::intersection(trianglePoly, clipMultiPoly, isect);

    // triangulate
    math::MultiPolygon isect2;
    isect2.reserve(isect.size());
    for (const auto &poly : isect) {
        isect2.push_back(ringPoints(poly.outer()));
        for (const auto &ring : poly.inners()) {
            isect2.push_back(ringPoints(ring));
        }
    }

    math::Triangles2d tris2(generalPolyTriangulate(isect2));

    // work around boost errorneously returning whole polygon as intersection
    // check area of input triangle against area of the result, should be same
    // or less. Definitelly it should not be significantly bigger.
    double interArea(0.0);
    for (const auto &t2 : tris2) {
        interArea += area(t2);
    }
    if (interArea > (area(tri2) * 1.1) ) { // 1.1 for numerical stability
        LOG(warn1) << "Throwing away spurious intersection (ratio of areas: "
                   << interArea / area(tri2) << ").";
        return {};
    }


    // restore Z coords
    math::Triangles3d tris3;
    tris3.reserve(tris2.size());
    for (const auto &t2 : tris2)
    {
        math::Triangle3d t3;
        for (int i = 0; i < 3; i++)
        {
            math::Point3 l(math::barycentricCoords(t2[i], tri2));
            t3[i](0) = t2[i](0);
            t3[i](1) = t2[i](1);
            t3[i](2) = l(0)*tri[0](2) +
                       l(1)*tri[1](2) +
                       l(2)*tri[2](2);
        }
        if (flip) {
            std::swap(t3[1], t3[2]);
        }
        tris3.push_back(t3);
    }

    return tris3;
}


math::Point3 barycentric3D(const math::Point3 &p, const math::Point3 &a,
                           const math::Point3 &b, const math::Point3 &c)
{
    typedef math::Point3 P3;
    double abp = norm_2(math::crossProduct(P3(b - a), P3(p - a)));
    double bcp = norm_2(math::crossProduct(P3(c - b), P3(p - b)));
    double cap = norm_2(math::crossProduct(P3(a - c), P3(p - c)));
    double nor = 1.0 / (abp + bcp + cap);
    return {bcp*nor, cap*nor, abp*nor};
}


std::tuple<math::Triangles3d, math::Triangles2d>
    clipTexturedTriangleNonconvex(const math::Triangle3d &tri,
                                  const math::Triangle2d &uv,
                                  const math::MultiPolygon &clipRegion)
{
    std::tuple<math::Triangles3d, math::Triangles2d> result;
    auto &tris3(std::get<0>(result));
    auto &uvs(std::get<1>(result));

    // clip the geometry first
    tris3 = clipTriangleNonconvex(tri, clipRegion);

    // if one same triangle -> copy texcoorsds
    // this clearly solves degen cases.
    if ((tris3.size() == 1) && (tris3[0] == tri)) {
        uvs.push_back(uv);
        return result;
    }

    // interpolate UV coords
    uvs.reserve(tris3.size());
    for (const auto &t3 : tris3)
    {
        math::Triangle2d t2;
        for (int i = 0; i < 3; i++)
        {
            math::Point3 l(barycentric3D(t3[i], tri[0], tri[1], tri[2]));
            for (int j = 0; j < 2; j++) {
                t2[i](j) = l(0)*uv[0](j) + l(1)*uv[1](j) + l(2)*uv[2](j);
            }
        }
        uvs.push_back(t2);
    }

    return result;
}

} // namespace geometry

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

//#include <boost/polygon/polygon.hpp>

#include "nonconvexclip.hpp"
#include "triangulate.hpp"

namespace geometry {

#if 0
namespace bp = boost::polygon;

typedef bp::point_data<int> Point;
typedef bp::polygon_with_holes_data<int> Polygon;
typedef bp::polygon_set_data<int> PolygonSet;

const double Multiplier = 1024;

math::Triangles3d clipTriangleNonconvex(const math::Triangle3d &tri,
                                        const Region &clipRegion)
{
    std::vector<Point> points;

    // convert the triangle to a Polygon
    points.reserve(3);
    for (const auto &p : tri) {
        points.emplace_back(p(0)*Multiplier, p(1)*Multiplier);
    }
    Polygon poly(points.begin(), points.end());

    // convert the clip region to a PolygonSet
    PolygonSet pset;
    for (const auto pts : clipRegion)
    {
        points.clear();
        points.reserve(pts.size());
        for (const auto& p : pts) {
            points.emplace_back(p(0)*Multiplier, p(1)*Multiplier);
        }
        pset.insert_vertex_sequence(points.begin(), points.end(),
                                    bp::COUNTERCLOCKWISE, false);
    }

    // calculate intersection
    std::vector<Polygon> clipped;
    PolygonSet intersection(pset & poly);
    intersection.get(clipped);
}
#endif

namespace bg = boost::geometry;

typedef bg::model::d2::point_xy<double> Point;
typedef bg::model::polygon<Point, false, false> Polygon;

template<typename List>
std::vector<Point> bgPoints(const List &list)
{
    std::vector<Point> result;
    result.reserve(list.size() + 1);
    for (const auto &p : list) {
        result.emplace_back(p(0), p(1));
    }
    return result;
}

math::Points2d outerPoints(const Polygon &poly)
{
    math::Points2d result;
    result.reserve(poly.outer().size());
    for (const auto &p : poly.outer()) {
        result.emplace_back(p.x(), p.y());
    }
    return result;
}

math::Triangles3d clipTriangleNonconvex(const math::Triangle3d &tri,
                                        const Region &clipRegion)
{
    // tri -> tri2
    math::Triangle2d tri2;
    for (int i = 0; i < 3; i++) {
        tri2[i](0) = tri[i](0);
        tri2[i](1) = tri[i](1);
    }

    // convert input to 2D polygons
    Polygon poly1, poly2;
    bg::append(poly1, bgPoints(tri));
    for (const auto &ring : clipRegion) {
        bg::append(poly2, bgPoints(ring));
    }

    // calculate intersection
    std::deque<Polygon> isect;
    bg::intersection(poly1, poly2, isect);

    // triangulate
    math::Triangles2d tris2;
    for (const auto &poly : isect) {
        auto tr(simplePolyTriangulate(outerPoints(poly)));
        tris2.insert(tris2.end(), tr.begin(), tr.end());
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
                                  const Region &clipRegion)
{
    std::tuple<math::Triangles3d, math::Triangles2d> result;
    auto &tris3(std::get<0>(result));
    auto &uvs(std::get<1>(result));

    // clip the geometry first
    tris3 = clipTriangleNonconvex(tri, clipRegion);

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

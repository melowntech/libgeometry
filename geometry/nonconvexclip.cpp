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

#include "nonconvexclip.hpp"
#include "triangulate.hpp"

namespace geometry {

namespace bg = boost::geometry;

typedef bg::model::d2::point_xy<double> Point;
typedef bg::model::polygon<Point, false, false> Polygon; // ccw, unclosed
typedef bg::model::multi_polygon<Polygon> MultiPolygon;
typedef Polygon::ring_type Ring;

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

math::Points2d ringPoints(const Ring &ring)
{
    math::Points2d result;
    result.reserve(ring.size());
    for (const auto &p : ring) {
        result.emplace_back(p.x(), p.y());
    }
    return result;
}

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
    double ccw = math::ccw(tri2[0], tri2[1], tri2[2]);

    if (std::abs(ccw) < 1e-10) {
        return {}; // TODO: handle exactly vertical triangles
    }

    // ensure counter-clockwise orientation for clipping
    if (ccw < 0.0) {
        std::swap(tri[1], tri[2]);
        std::swap(tri2[1], tri2[2]);
        flip = true;
    }

    // convert input to 2D polygons
    Polygon poly1;
    bg::assign_points(poly1, bgPoints(tri));

    MultiPolygon poly2;
    for (const auto &pts : clipRegion) {
        Polygon part;
        bg::assign_points(part, bgPoints(pts));
        poly2.push_back(part);
    }

    // calculate intersection
    std::deque<Polygon> isect;
    bg::intersection(poly1, poly2, isect);

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

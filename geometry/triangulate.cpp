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

#include "boost/polygon/voronoi.hpp"

#include "dbglog/dbglog.hpp"

#include "triangulate.hpp"

namespace geometry {

namespace {

const double eps = 1e-12;

typedef math::Point2d Point;

double convex(const Point &a, const Point &b, const Point &c)
{
    return (b(0) - a(0))*(c(1) - a(1)) - (b(1) - a(1))*(c(0) - a(0));
}

bool insideTriangle(const Point &a, const Point &b, const Point &c,
                    const Point &pt)
{
    return convex(a, b, pt) > eps &&
           convex(b, c, pt) > eps &&
           convex(c, a, pt) > eps;
}

struct Vertex
{
    Point pt;
    double cx;
    bool ear;

    Vertex(const Point &pt) : pt(pt), cx(0), ear(false) {}
};

typedef std::vector<Vertex> Vertices;

bool isEar(const Vertices &vert, int i, int j, int k)
{
    // the vertex must be convex
    if (convex(vert[i].pt, vert[j].pt, vert[k].pt) < eps) {
        return false;
    }

    // the remaining points must not lie inside the i-j-k triangle
    for (int l = 0; l < int(vert.size()); l++)
    {
        if (l != i && l != j && l != k)
        {
            if (insideTriangle(vert[i].pt, vert[j].pt, vert[k].pt, vert[l].pt)) {
                return false;
            }
        }
    }
    return true;
}

} // namespace

math::Triangles2d simplePolyTriangulate(const math::Polygon &polygon)
{
    std::vector<Vertex> vert;
    vert.reserve(polygon.size());
    for (const auto &pt : polygon) {
        vert.push_back(pt);
    }

    math::Triangles2d result;

    int j = 0, loop = 0;
    while (vert.size() >= 3)
    {
        if (j >= int(vert.size())) {
            j = 0;
        }
        int n = vert.size();
        int i = j-1, k = j+1;
        if (i < 0) { i = n-1; }
        if (k >= n) { k = 0; }

        // TODO: process ears with shorter diagnoals first, for better quality

        if (isEar(vert, i, j, k))
        {
            result.push_back({{vert[i].pt, vert[j].pt, vert[k].pt}});
            vert.erase(vert.begin() + j);
            loop = 0;
        }
        else
        {
            j++;
            if (++loop >= n) {
                LOG(err2) << "Can't triangulate polygon, not simple or not ccw?";
                break;
            }
        }
    }
    return result;
}

namespace {

// based on https://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon
// TODO: maybe move to polygon.hpp
bool pointInPolygon(const math::Point2d &test, const math::Polygon &poly)
{
    bool c = false;
    int nvert = poly.size();
    for (int i = 0, j = nvert-1; i < nvert; j = i++)
    {
        if ( ((poly[i](1) > test(1)) != (poly[j](1) > test(1))) &&
             (test(0) < (poly[j](0) - poly[i](0)) * (test(1) - poly[i](1)) /
                        (poly[j](1) - poly[i](1)) + poly[i](0)) )
        {
           c = !c;
        }
    }
    return c;
}

bool pointInMultiPolygon(const math::Point2d &test,
                         const math::MultiPolygon &mpoly)
{
    bool c = false;
    for (const auto &poly : mpoly)
    {
        int nvert = poly.size();
        for (int i = 0, j = nvert-1; i < nvert; j = i++)
        {
            if ( ((poly[i](1) > test(1)) != (poly[j](1) > test(1))) &&
                 (test(0) < (poly[j](0) - poly[i](0)) * (test(1) - poly[i](1)) /
                            (poly[j](1) - poly[i](1)) + poly[i](0)) )
            {
               c = !c;
            }
        }
    }
    return c;
}

} // namespace

namespace bp = boost::polygon;

const double Multiplier = (1 << 16);
const double InvMultiplier = 1.0 / Multiplier;


math::Triangles2d generalPolyTriangulate(const math::MultiPolygon &mpolygon)
{
    std::vector<bp::point_data<int> > vertices;
    for (const auto &poly : mpolygon) {
        for (const auto &p : poly) {
            vertices.emplace_back(p(0)*Multiplier, p(1)*Multiplier);
        }
    }

    bp::voronoi_diagram<double> voronoi;
    bp::construct_voronoi(vertices.begin(), vertices.end(), &voronoi);

    math::Triangles2d result;

    for (const auto& vertex : voronoi.vertices())
    {
        std::vector<math::Point2d> tri;
        auto edge = vertex.incident_edge();
        do
        {
            auto cell = edge->cell();
            assert(cell->contains_point());

            const auto &v = vertices[cell->source_index()];
            tri.emplace_back(v.x()*InvMultiplier, v.y()*InvMultiplier);

            if (tri.size() == 3)
            {
                math::Point2d c = (tri[0] + tri[1] + tri[2]) * (1.0 / 3);
                if (pointInMultiPolygon(c, mpolygon))
                {
                    result.push_back({{tri[0], tri[1], tri[2]}});
                }
                tri.erase(tri.begin() + 1);
            }

            edge = edge->rot_next();
        }
        while (edge != vertex.incident_edge());
    }
    return result;
}

} // namespace geometry

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

#include "utility/gccversion.hpp"

#include "polygon.hpp"

#include "triangulate.hpp"

#include "CDT/CDT.h"

namespace geometry {

namespace bp = boost::polygon;

namespace {

const double eps = 1e-12;

typedef math::Point2d Point;

bool insideTriangle(const Point &a, const Point &b, const Point &c,
                    const Point &pt)
{
    return math::ccw(a, b, pt) > eps &&
           math::ccw(b, c, pt) > eps &&
           math::ccw(c, a, pt) > eps;
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
    if (math::ccw(vert[i].pt, vert[j].pt, vert[k].pt) < eps) {
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
        auto n = int(vert.size());
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

#if 0
// based on https://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon
// TODO: maybe move to polygon.hpp
bool pointInPolygon(const math::Point2d &test, const math::Polygon &poly)
    UTILITY_POSSIBLY_UNUSED;
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
#endif

// Converts input coordinates to interval [0, 2^16] and vice versa
// (work around boost polygon integer nature)
class CoordConverter {
    math::Point2 offset_;
    double mult_;
    double invMult_;

public:
    explicit CoordConverter(const math::MultiPolygon& mpoly)
    {
        math::Extents2 ext(math::InvalidExtents{});
        for (const auto& poly : mpoly) {
            for (const auto& p : poly) {
                math::update(ext, p);
            }
        }
        auto size = math::size(ext);
        const double MaxValue = 1 << 16;

        mult_ = MaxValue / std::max(size.width, size.height);
        invMult_ = 1.0 / mult_;
        offset_ = ext.ll;
    }

    inline bp::point_data<double> operator()(const math::Point2& p) const
    {
        return bp::point_data<double>((p(0) - offset_(0)) * mult_,
                                      (p(1) - offset_(1)) * mult_);
    }

    inline math::Point2 operator()(const bp::point_data<double>& p) const
    {
        return math::Point2(p.x() * invMult_ + offset_(0),
                            p.y() * invMult_ + offset_(1));
    }
};

} // namespace

math::Triangles2d generalPolyTriangulate(const math::MultiPolygon &mpolygon)
{
    CoordConverter cnv(mpolygon);
    std::vector<bp::point_data<double> > points;
    for (const auto &poly : mpolygon) {
        for (const auto &p : poly) {
            points.emplace_back(cnv(p));
        }
    }

    bp::voronoi_diagram<double> voronoi;
    bp::construct_voronoi(points.begin(), points.end(), &voronoi);

    math::Triangles2d result;

    // TODO: we should really do constrained Delaunay triangulation here to
    // make sure all required edges are contained in the output.

    for (const auto& vertex : voronoi.vertices())
    {
        std::vector<math::Point2d> tri;
        auto edge = vertex.incident_edge();
        do
        {
            auto cell = edge->cell();
            assert(cell->contains_point());

            const auto &v = points[cell->source_index()];
            tri.emplace_back(cnv(v));

            if (tri.size() == 3)
            {
                math::Point2d c = (tri[0] + tri[1] + tri[2]) * (1.0 / 3);
                if (insideMultiPolygon(mpolygon, c))
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

math::Triangles2d generalPolyTriangulateCDT(const math::MultiPolygon &mpolygon)
{
    std::vector<CDT::V2d<double>> cdtVertices;
    std::vector<CDT::Edge> cdtEdges;

    // Note math::MultiPolygon consists of open CCW rings and open CW holes
    // (ie. first point of the ring is not same as the last point)

    for (const auto &poly : mpolygon) {
        if (poly.size()) {
            CDT::VertInd cdtVertexOffset = static_cast<CDT::VertInd>(cdtVertices.size());
            CDT::VertInd cdtVertexId = 0;
            cdtVertices.reserve(cdtVertices.size() + poly.size());
            cdtEdges.reserve(cdtEdges.size() + poly.size() - 1);
            for (const auto &p : poly) {
                cdtVertices.push_back({p(0), p(1)});
                cdtEdges.emplace_back(
                    cdtVertexOffset + cdtVertexId, 
                    cdtVertexOffset + (cdtVertexId + 1) % static_cast<CDT::VertInd>(poly.size()));
                cdtVertexId++;
            }
        }
    }

    if (cdtVertices.size() == 0 || cdtEdges.size() == 0) {
        return {};
    }

    [[maybe_unused]] CDT::DuplicatesInfo di =
        CDT::RemoveDuplicatesAndRemapEdges(cdtVertices, cdtEdges);

    // Pre-conditions:
    // - No duplicated points (use provided functions for removing duplicate points and re-mapping edges)
    // - No two constraint edges intersect each other (overlapping boundaries are allowed)
    // Post-conditions:
    // - Triangles have counter-clockwise (CCW) winding
    
    CDT::Triangulation cdt(
        CDT::VertexInsertionOrder::Auto,
        CDT::IntersectingConstraintEdges::Resolve,
        /*minDistToConstraintEdge*/ 0.0);

    cdt.insertVertices(cdtVertices);
    cdt.insertEdges(cdtEdges);
    cdt.eraseOuterTrianglesAndHoles();

    math::Triangles2d tris2d;
    tris2d.reserve(cdt.triangles.size());
    for (const auto & cdtTri : cdt.triangles) {
        tris2d.emplace_back(math::Triangle2d{
            math::Point2d{cdt.vertices[cdtTri.vertices[0]].x, cdt.vertices[cdtTri.vertices[0]].y},
            math::Point2d{cdt.vertices[cdtTri.vertices[1]].x, cdt.vertices[cdtTri.vertices[1]].y},
            math::Point2d{cdt.vertices[cdtTri.vertices[2]].x, cdt.vertices[cdtTri.vertices[2]].y}
        });
    }
    return tris2d;
}

} // namespace geometry

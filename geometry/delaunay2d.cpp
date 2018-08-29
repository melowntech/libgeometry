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

// This is a workaround for CGAL Assertions currently ruining debug mode because
// CGAL expects flags which are not currently set
#ifndef NDEBUG
#  define NDEBUG
#endif

// WARNING: CGAL is GPL
// TODO: consider using the GNU Triangulated Surface Library, which is LGPL

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include "delaunay2d.hpp"

namespace geometry {

std::vector<DTriangle> delaunayTriangulation2d(const math::Points2 &points)
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb>                 Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds>                   Delaunay;
    typedef K::Point_2                                               Point;

    std::vector<std::pair<Point, unsigned> > cgalPoints;
    {
        unsigned index = 0;
        for (const auto &pt : points) {
            cgalPoints.emplace_back(Point(pt(0), pt(1)), index++);
        }
    }

    Delaunay triangulation;
    triangulation.insert(cgalPoints.begin(), cgalPoints.end());

    std::vector<DTriangle> triangles;
    triangles.reserve(triangulation.number_of_faces());

    for (auto fit = triangulation.finite_faces_begin();
              fit != triangulation.finite_faces_end(); ++fit)
    {
        DTriangle indices;
        for (int i = 0; i < 3; i++) {
            indices[i] = fit->vertex(i)->info();
        }
        triangles.push_back(indices);
    }

    return triangles;
}

#ifdef GEOMETRY_HAS_CGAL_4_11
void constrainedDelaunayTriangulation2d(
        const math::Points2 &points,
        const std::vector<DEdge> &constrained_edges,
        math::Points2 &out_points,
        std::vector<DTriangle> &triangles)
{
    // adapted from
    // cgal/Triangulation_2/examples/Triangulation_2/polylines_triangulation.cpp

    typedef CGAL::Exact_predicates_exact_constructions_kernel                  K;
    typedef CGAL::Exact_intersections_tag                                   Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K,CGAL::Default,Itag> CDT;
    typedef CGAL::Constrained_triangulation_plus_2<CDT>                     CDTP;

    std::vector<K::Point_2> cgalPoints;
    cgalPoints.reserve(points.size());
    for (const auto &p : points) {
        cgalPoints.emplace_back(p(0), p(1));
    }

    // build triangulation
    CDTP cdtp;
    for (const auto &ce : constrained_edges)
    {
        cdtp.insert_constraint(cgalPoints[ce[0]], cgalPoints[ce[1]]);
    }
    cdtp.insert(cgalPoints.begin(), cgalPoints.end());

    // return new points
    out_points.clear();
    out_points.reserve(cdtp.tds().vertices().size());

    for (const auto &v : cdtp.tds().vertices())
    {
        out_points.emplace_back(CGAL::to_double(v.point().x()),
                                CGAL::to_double(v.point().y()));
    }

    // return triangles
    triangles.clear();
    triangles.reserve(cdtp.number_of_faces());

    for (auto fit = cdtp.finite_faces_begin();
              fit != cdtp.finite_faces_end(); ++fit)
    {
        DTriangle indices;
        for (int i = 0; i < 3; i++) {
            indices[i] = cdtp.tds().vertices().index(fit->vertex(i));
        }
        triangles.push_back(indices);
    }
}
#endif

} // namespace geometry

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
#include "dbglog/dbglog.hpp"
#include "gtsutils.hpp"
#include "utility/expect.hpp"
#include "utility/small_map.hpp"

#include "geometry/polygon.hpp"

#include "multipolygon-triangulate-gts.hpp"

namespace geometry
{
/// Helper struct for triangle extraction
struct TriData
{
    const utility::small_map<GtsVertex*, ItPair>& vertexMap;
    const math::MultiPolygon& mpolygon;
    std::vector<TriangleItPair> triangles;
};

/// Helper to find vertex in vertexMap and throw error if not found
inline ItPair
    findVertexInMap(GtsVertex* v,
                    const utility::small_map<GtsVertex*, ItPair>& vertexMap)
{
    auto it = vertexMap.find(v);
    if (it == vertexMap.end())
    {
        LOGTHROW(err4, std::runtime_error)
            << "Cannot find a vertex in added vertex map.";
    }
    return it->second;
}

/// Extracts triangle and adds it to `data.triangles` if center is in
/// multipolygon
gint extractTriangle(GtsTriangle* t, TriData* data)
{
    GtsVertex *v1, *v2, *v3;
    gts_triangle_vertices(t, &v1, &v2, &v3);

    // check that centroid is inside multipolygon
    math::Point2 p1(v1->p.x, v1->p.y);
    math::Point2 p2(v2->p.x, v2->p.y);
    math::Point2 p3(v3->p.x, v3->p.y);
    math::Point2d c = (p1 + p2 + p3) * (1.0 / 3);

    if (geometry::insideMultiPolygon(data->mpolygon, c))
    {
        data->triangles.push_back({ findVertexInMap(v1, data->vertexMap),
                                    findVertexInMap(v2, data->vertexMap),
                                    findVertexInMap(v3, data->vertexMap) });
    }
    return {};
}

std::vector<TriangleItPair>
    multipolygonTriangulateGts(const math::MultiPolygon& mpolygon)
{
    geometry::checkGtsInitialized();

    GtsSurface* s = gts_surface_new(gts_surface_class(),
                                    gts_face_class(),
                                    gts_edge_class(),
                                    gts_vertex_class());

    // create vertices
    GSList* verticesList = nullptr;
    utility::small_map<GtsVertex*, ItPair> vertexMap;
    for (auto polyIt = mpolygon.begin(); polyIt != mpolygon.end(); ++polyIt)
    {
        for (auto ptIt = polyIt->begin(); ptIt != polyIt->end(); ++ptIt)
        {
            GtsVertex* vx = gts_vertex_new(gts_vertex_class(),
                                           (*ptIt)(0),
                                           (*ptIt)(1),
                                           0.0);
            verticesList = g_slist_prepend(verticesList, vx);
            vertexMap[vx] = std::make_pair(polyIt, ptIt);
        }
    }
    // prepend & reverse is more efficient for linked-lists
    verticesList = g_slist_reverse(verticesList);

    // add boundary face enclosing all vertices
    GtsTriangle* t
        = gts_triangle_enclosing(gts_triangle_class(), verticesList, 1.1);
    GtsVertex *tv1, *tv2, *tv3;
    gts_triangle_vertices(t, &tv1, &tv2, &tv3);
    gts_surface_add_face(s,
                         gts_face_new(gts_face_class(), t->e1, t->e2, t->e3));

    // add vertices to surface
    GSList* l = verticesList;
    while (l)
    {
        GtsVertex* ret = gts_delaunay_add_vertex(s, (GtsVertex*)l->data, NULL);
        if(ret == l->data)
        {
            LOGTHROW(err4, std::runtime_error)
                << "Vertex to be added is not contained in the convex hull "
                   "bounding surface, but it should be (triangle enclosing all "
                   "vertices was added).";
        }
        if(ret != nullptr)
        {
            LOG(warn2) << "Found duplicated vertex in the polygon (can be "
                          "holes touching)";
            // replace with the first vertex
            gts_object_destroy(GTS_OBJECT(l->data));
            l->data = ret;
        }

        l = l->next;
    }

    // add constraints for boundaries of all polygons
    l = verticesList;
    for (auto& poly : mpolygon)
    {
        GtsVertex* firstVertex = (GtsVertex*)l->data; // save first vertex

        // add constraints for all but last vertex
        for (std::size_t i = 0; i < (poly.size() - 1); ++i)
        {
            utility::expect(l != nullptr, "Unexpected end of vertex list");
            utility::expect(l->next != nullptr,
                            "Unexpected end of vertex list");

            GtsConstraint* c = (GtsConstraint*)gts_edge_new(
                GTS_EDGE_CLASS(gts_constraint_class()),
                (GtsVertex*)l->data,
                (GtsVertex*)l->next->data);

            utility::expect(gts_delaunay_add_constraint(s, c) == NULL,
                            "Unable to add constraint");

            l = l->next;
        }
        utility::expect(l != nullptr, "Unexpected end of vertex list");

        // add constraint last-first
        GtsConstraint* c = (GtsConstraint*)gts_edge_new(
            GTS_EDGE_CLASS(gts_constraint_class()),
            (GtsVertex*)l->data,
            firstVertex);
        utility::expect(gts_delaunay_add_constraint(s, c) == NULL,
                        "Unable to add constraint");

        l = l->next; // vertex of next polygon
    }
    utility::expect(
        l == nullptr,
        "Vertex list is longer than number of vertices in mpolygon");

    // remove the boundary face
    gts_allow_floating_vertices = TRUE;
    gts_object_destroy(GTS_OBJECT(tv1));
    gts_object_destroy(GTS_OBJECT(tv2));
    gts_object_destroy(GTS_OBJECT(tv3));
    gts_allow_floating_vertices = FALSE;

    // get the faces
    TriData tridata { vertexMap, mpolygon, {} }; // helper data
    gts_surface_foreach_face(s, (GtsFunc)extractTriangle, &tridata);

    // cleanup
    gts_object_destroy(GTS_OBJECT(s));

    return tridata.triangles;
}

} // namespace geometry

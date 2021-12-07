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

#include "delaunay2d.hpp"
#include "gtsutils.hpp"
#include "meshop.hpp"

namespace geometry {

namespace {

const unsigned SUPER_IDX = unsigned(-1);

GtsTriangle* enclosingTriangle(const math::Points2& points)
{
    GSList* list = nullptr;
    for (const math::Point2& p : points) {
        GtsPoint* point = gts_point_new(gts_point_class(), p(0), p(1), 0);
        list = g_slist_prepend(list, point);
    }
    GtsTriangle* t = gts_triangle_enclosing(gts_triangle_class(), list, 10.);
    GtsVertex *v1, *v2, *v3;
    gts_triangle_vertices(t, &v1, &v2, &v3);
    GTS_OBJECT(v1)->reserved = GINT_TO_POINTER(SUPER_IDX);
    GTS_OBJECT(v2)->reserved = GINT_TO_POINTER(SUPER_IDX);
    GTS_OBJECT(v3)->reserved = GINT_TO_POINTER(SUPER_IDX);
    g_slist_free(list);
    return t;
}

gint extractFace(GtsTriangle *t, std::vector<DTriangle>* triangles)
{
    GtsVertex *v1, *v2, *v3;
    gts_triangle_vertices(t, &v1, &v2, &v3);
    // retrieve the index from reserved
    const unsigned i1 = GPOINTER_TO_UINT(GTS_OBJECT(v1)->reserved);
    const unsigned i2 = GPOINTER_TO_UINT(GTS_OBJECT(v2)->reserved);
    const unsigned i3 = GPOINTER_TO_UINT(GTS_OBJECT(v3)->reserved);
    if (i1 != SUPER_IDX && i2 != SUPER_IDX && i3 != SUPER_IDX) {
        triangles->emplace_back(DTriangle{ i1, i2, i3 });
    }
    return 0;
}

} // namespace

std::vector<DTriangle> delaunayTriangulation2d(const math::Points2 &points)
{
    checkGtsInitialized();
    GtsSurface* gts = gts_surface_new(gts_surface_class(),
                                      gts_face_class(),
                                      gts_edge_class(),
                                      gts_vertex_class());
    GtsTriangle* super = enclosingTriangle(points);
    GtsFace* f =
        gts_face_new(gts_face_class(), super->e1, super->e2, super->e3);
    gts_surface_add_face(gts, f);

    for (std::size_t i = 0; i < points.size(); ++i) {
        const math::Point2& p = points[i];
        GtsVertex* v = gts_vertex_new(gts_vertex_class(), p(0), p(1), 0);
        // store the point index to reserved
        GTS_OBJECT(v)->reserved = GUINT_TO_POINTER(i);
        GtsVertex* v1 = gts_delaunay_add_vertex(gts, v, nullptr);
        g_assert(v1 != v);
    }

    std::vector<DTriangle> triangles;
    gts_surface_foreach_face(gts, (GtsFunc)extractFace, &triangles);

    gts_object_destroy(GTS_OBJECT(gts));
    return triangles;
}

} // namespace geometry

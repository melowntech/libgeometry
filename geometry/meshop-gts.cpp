/**
 * Copyright (c) 2018 Melown Technologies SE
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
 * @file meshop-gts.cpp
 * @author Alena Chernikava  <alena.chernikava@melown.com>
 *
 * simplification with Lindstrom-Turk algorithm from GTS library
 */

#include <gts.h>
#include "./meshop.hpp"
#include "utility/openmp.hpp"

namespace geometry {

// Looks for an existing gts edge.
// if it is found then return it
// If it is not found -> insert into a map of edges and return it
// Just a helper
GtsEdge * get_gts_edge (std::map <std::pair<int, int>, GtsEdge *> &edges
    , const std::vector <GtsVertex *> &vertices
    , int v1
    , int v2)
{
    int vmin = std::min (v1, v2);
    int vmax = std::max (v1, v2);
    auto e_lookup = edges.find(std::make_pair (vmin, vmax));
    if (e_lookup == edges.end()) {
        GtsEdge *e = gts_edge_new (gts_edge_class (), vertices[vmin], vertices[vmax]);
        edges.emplace (std::make_pair (vmin, vmax), e);
        return e;
    }
    return e_lookup->second;
}

// Convert vadstena::Mesh to gts::Mesh
void toGTSMesh(const geometry::Mesh &mesh, GtsSurface &s)
{
    std::vector <GtsVertex*> gts_vertices;
    gts_vertices.reserve(mesh.vertices.size());
    for (const auto& v : mesh.vertices) {
        GtsVertex *vx = gts_vertex_new (gts_vertex_class (), v(0), v(1), v(2));
        gts_vertices.push_back(vx);
    }

    std::map <std::pair<int, int>, GtsEdge *> gts_edges;
    // create GTS faces
    for (const geometry::Face& face : mesh.faces)
    {
        // Build edges
        GtsEdge *e1 = get_gts_edge (gts_edges, gts_vertices, face.a, face.b);
        GtsEdge *e2 = get_gts_edge (gts_edges, gts_vertices, face.b, face.c);
        GtsEdge *e3 = get_gts_edge (gts_edges, gts_vertices, face.c, face.a);
        // Build the face
        GtsFace *f = gts_face_new (gts_face_class (), e1, e2, e3);
        gts_surface_add_face (&s, f);
    }
}

// A 'callback' for foreach_vertex call that converts
// all verticies from gts to vadstena and add it to the list of vertices
// Just a helper
void vertexGts2Va (GtsPoint *gts_p, gpointer *data)
{
    math::Points3d *va_points = (math::Points3d *) data[0];
    uint *current_index = (uint *) data[1];
    va_points->emplace_back(gts_p->x, gts_p->y, gts_p->z);
    // Here we save index of vertex
    // This is inspired by gts codebase
    GTS_OBJECT (gts_p)->reserved = GUINT_TO_POINTER (*current_index);
    ++(*current_index);
}

// A 'callback' for foreach_face call that converts
// all faces from gts to vadstena and add it to the list of faces
// Just a helper
void faceGts2Va (GtsTriangle *t, Face::list *va_faces)
{
  GtsVertex *v1, *v2, *v3;
  gts_triangle_vertices (t, &v1, &v2, &v3);
  // Create our face with indexes of verteces!
  va_faces->emplace_back(
          GPOINTER_TO_UINT (GTS_OBJECT (v1)->reserved),
          GPOINTER_TO_UINT (GTS_OBJECT (v2)->reserved),
          GPOINTER_TO_UINT (GTS_OBJECT (v3)->reserved));
}

// Convert gts::Mesh to vadstena::Mesh
void fromGts(geometry::Mesh &mesh, GtsSurface &s)
{
    // create our vertices
    uint vertex_index = 0;
    gpointer data[2];
    data[0] = &mesh.vertices;
    data[1] = &vertex_index;
    gts_surface_foreach_vertex (&s, (GtsFunc) vertexGts2Va, data);
    gts_surface_foreach_face (&s, (GtsFunc) faceGts2Va, &mesh.faces);
}

// TODO faceCount
geometry::Mesh::pointer simplify_gts(const geometry::Mesh &mesh, long edgeCountMax)
{
    // the make_gts_class_system_threadsafe (void);
    // MUST be called in main thread! before this moment
    GtsSurface *s = gts_surface_new (gts_surface_class (),
            gts_face_class (),
            gts_edge_class (),
            gts_vertex_class ());
    toGTSMesh(mesh, *s);

    //  set maximum fold angle to F degrees, default is one degree"
    gdouble fold = M_PI/180.;
    //  fold = atof (optarg)*PI/180.;
    GtsStopFunc stop_func =  (GtsStopFunc) gts_coarsen_stop_number;
    GtsCoarsenFunc coarsen_func = (GtsCoarsenFunc) gts_volume_optimized_vertex;
    // Default parameters
    GtsVolumeOptimizedParams params = { 0.5, 0.5, 0. };
    // TODO export to Window-fusion
    //params.volume_weight = atof (optarg);
    //params.boundary_weight = atof (optarg);
    //params.shape_weight = atof (optarg);

    GtsKeyFunc cost_func = (GtsKeyFunc) gts_volume_optimized_cost;
    gpointer coarsen_data = &params;
    gpointer cost_data = &params;

    gts_surface_coarsen (s,
			 cost_func, cost_data,
			 coarsen_func, coarsen_data,
			 stop_func, &edgeCountMax, fold);

    auto newMesh(std::make_shared<geometry::Mesh>());
    fromGts (*newMesh, *s);
    gts_object_destroy (GTS_OBJECT (s));
    return newMesh;
}

void getSubSurface (void *ee, gpointer *data)
{
    GtsSurface *s_sub = (GtsSurface *) data[1];

    math::Extents2 *cellExtents = (math::Extents2 *) data[0];

    GtsFace *f = (GtsFace *) ee;
    GtsVertex *v1, *v2, *v3;
    gts_triangle_vertices (GTS_TRIANGLE (f), &v1, &v2, &v3);

    gint result = 1;

    gdouble midx = GTS_POINT (v1) -> x + GTS_POINT (v2)->x + GTS_POINT (v3) ->x;
    midx /= 3;
    gdouble midy = GTS_POINT (v1) -> y + GTS_POINT (v2)->y + GTS_POINT (v3) ->y;
    midy /= 3;
    if (midx < cellExtents->ll(0)) {
        result = 0;
    } else if (midx > cellExtents->ur(0)) {
        result = 0;
    } else if (midy < cellExtents->ll(1)) {
        result = 0;
    } else if (midy > cellExtents->ur(1)) {
        result = 0;
    }
    if (result) {
        gts_surface_add_face (s_sub, f);
    }
}

void make_gts_class_system_threadsafe (void) {
    // Claass system in GTS is not thread safe at compile time
    // initialise it in the main thread -> then it would thread-safe at run time
    gts_face_class ();
    gts_edge_class ();
    gts_bbox_class ();
    gts_graph_class ();
    gts_split_class ();
    gts_point_class ();
    gts_gnode_class ();
    gts_pnode_class ();
    gts_fnode_class ();
    gts_gedge_class ();
    gts_nedge_class ();
    gts_nface_class ();
    gts_pgedge_class ();
    gts_wgedge_class ();
    gts_object_class ();
    gts_hsplit_class ();
    gts_vertex_class ();
    gts_pgraph_class ();
    gts_ngnode_class ();
    gts_wgnode_class ();
    gts_wgraph_class ();
    gts_cluster_class ();
    gts_surface_class ();
    gts_nvertex_class ();
    gts_segment_class ();
    gts_hsurface_class ();
    gts_triangle_class ();
    gts_psurface_class ();
    gts_list_face_class ();
    gts_containee_class ();
    gts_container_class ();
    gts_constraint_class ();
    gts_gnode_split_class ();
    gts_color_vertex_class ();
    gts_cluster_grid_class ();
    gts_vertex_normal_class ();
    gts_surface_inter_class ();
    gts_hash_container_class ();
    gts_slist_containee_class ();
    gts_slist_container_class ();
}

geometry::Mesh::pointer simplify_gts_in_grid(const geometry::Mesh &mesh
    , std::vector<std::vector <geometry::GridCell>> &gridCells, bool inParallel)
{
    // the make_gts_class_system_threadsafe (void);
    // MUST be called in main thread! be fore this moment
    GtsSurface *s = gts_surface_new (gts_surface_class (),
            gts_face_class (),
            gts_edge_class (),
            gts_vertex_class ());
    toGTSMesh(mesh, *s);

    //  set maximum fold angle to F degrees, default is one degree"
    gdouble fold = M_PI/180.;
    //  fold = atof (optarg)*PI/180.;
    GtsStopFunc stop_func =  (GtsStopFunc) gts_coarsen_stop_number;
    GtsCoarsenFunc coarsen_func = (GtsCoarsenFunc) gts_volume_optimized_vertex;
    // Default parameters
    // TODO
    GtsVolumeOptimizedParams params = { 0.5, 0.5, 0. };
    //params.volume_weight = atof (optarg);
    //params.boundary_weight = atof (optarg);
    //params.shape_weight = atof (optarg);

    GtsKeyFunc cost_func = (GtsKeyFunc) gts_volume_optimized_cost;
    gpointer coarsen_data = &params;
    gpointer cost_data = &params;

    std::map <std::pair <long, long>, GtsSurface*> subSurfaces;
    auto rows = gridCells.size();
    auto columns = gridCells[0].size(); // All rows have the same number of columns
    auto oneCellProcess ([&](long row, long column, GtsSurface *subSurface) -> void
    {
        auto edgeMaxCount = gridCells[row][column].maxFaceCount * 3 / 2;
        LOG (info3) << "Simplifying cell [" << row << "," << column << "]" << "/[" << rows << "," << columns << "] "
            << gridCells[row][column].extent
            << " to " << gridCells[row][column].maxFaceCount << " faces";
        gts_surface_coarsen (subSurface,
                cost_func, cost_data,
                coarsen_func, coarsen_data,
                stop_func, &edgeMaxCount, fold);

        gts_object_destroy (GTS_OBJECT (subSurface));
    });


    if (!inParallel) {
        for (uint row = 0 ; row < rows; row+=1)
        for (uint column = 0 ; column < columns; column+=1)
        {
            GtsSurface *subSurface = gts_surface_new (gts_surface_class (),
                    gts_face_class (),
                    gts_edge_class (),
                    gts_vertex_class ());
            gpointer subParams[2];
            subParams[0] = &gridCells[row][column].extent;
            subParams[1] = subSurface;

            // Fill the sub surface with faces from this cell and calculate mesh_area in cells
            gts_surface_foreach_face (s, (GtsFunc) getSubSurface, &subParams);
            oneCellProcess (row, column, subSurface);
        }
    } else {
        auto computeBunchInParallel ([&](long cellRowStart, long cellColumnStart, long cellRowShift, long cellColumnShift) -> void
        {
            // We cannot change a gts_surface (simplifying one cell of the mesh)
            // and at the same time iterate over gts_surface_verteces/faces/edges
            // (creating a subsurface for another cell)
            // The error:
            //     Gts:ERROR:surface.c:184:gts_surface_remove_face: assertion failed: (s->keep_faces == FALSE)
            // Solution: prepare all sub surfaces for this bunch in advance

            for (uint row = cellRowStart ; row < rows; row+=cellRowShift)
            for (uint column = cellColumnStart ; column < columns; column+=cellColumnShift)
            {
                GtsSurface *subSurface = gts_surface_new (gts_surface_class (),
                        gts_face_class (),
                        gts_edge_class (),
                        gts_vertex_class ());
                gpointer subParams[2];
                subParams[0] = &gridCells[row][column].extent;
                subParams[1] = subSurface;

                // Fill the sub surface with faces from this cell
                gts_surface_foreach_face (s, (GtsFunc) getSubSurface, &subParams);
                subSurfaces.emplace (std::make_pair <long, long> (row, column), subSurface);
            }

            UTILITY_OMP(parallel for schedule(dynamic, 1) collapse(2))
            for (uint row = cellRowStart ; row < rows; row+=cellRowShift)
            for (uint column = cellColumnStart ; column < columns; column+=cellColumnShift)
            {
                GtsSurface *subSurface = subSurfaces.find (std::make_pair <long, long> (row, column))->second;
                oneCellProcess (row, column, subSurface);
            }
        });
        computeBunchInParallel (0, 0, 2, 2);
        subSurfaces.clear();
        computeBunchInParallel (0, 1, 2, 2);
        subSurfaces.clear();
        computeBunchInParallel (1, 0, 2, 2);
        subSurfaces.clear();
        computeBunchInParallel (1, 1, 2, 2);
    }

    auto newMesh(std::make_shared<geometry::Mesh>());
    fromGts (*newMesh, *s);
    gts_object_destroy (GTS_OBJECT (s));
    return newMesh;
}

} // namespace end

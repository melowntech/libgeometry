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

#include "gtsutils.hpp"
#include "meshop.hpp"
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
gint vertexGts2Va (GtsPoint *gts_p, gpointer *data)
{
    math::Points3d *va_points = (math::Points3d *) data[0];
    unsigned int *current_index = (unsigned int *) data[1];
    va_points->emplace_back(gts_p->x, gts_p->y, gts_p->z);
    // Here we save index of vertex
    // This is inspired by gts codebase
    GTS_OBJECT (gts_p)->reserved = GUINT_TO_POINTER (*current_index);
    ++(*current_index);

    return {};
}

// A 'callback' for foreach_face call that converts
// all faces from gts to vadstena and add it to the list of faces
// Just a helper
gint faceGts2Va (GtsTriangle *t, Face::list *va_faces)
{
  GtsVertex *v1, *v2, *v3;
  gts_triangle_vertices (t, &v1, &v2, &v3);
  // Create our face with indexes of verteces!
  va_faces->emplace_back(
          GPOINTER_TO_UINT (GTS_OBJECT (v1)->reserved),
          GPOINTER_TO_UINT (GTS_OBJECT (v2)->reserved),
          GPOINTER_TO_UINT (GTS_OBJECT (v3)->reserved));

    return {};
}

// Convert gts::Mesh to vadstena::Mesh
void fromGts(geometry::Mesh &mesh, GtsSurface &s)
{
    // create our vertices
    unsigned int vertex_index = 0;
    gpointer data[2];
    data[0] = &mesh.vertices;
    data[1] = &vertex_index;
    gts_surface_foreach_vertex (&s, (GtsFunc) vertexGts2Va, data);
    gts_surface_foreach_face (&s, (GtsFunc) faceGts2Va, &mesh.faces);
}

// TODO faceCount
geometry::Mesh::pointer simplify_gts(const geometry::Mesh &mesh
                                     , long long edgeCountMax, double costMax)
{
    // GTS class system MUST be initialized from the main thread before this moment
    checkGtsInitialized();
    GtsSurface *s = gts_surface_new (gts_surface_class (),
            gts_face_class (),
            gts_edge_class (),
            gts_vertex_class ());
    toGTSMesh(mesh, *s);

    //  set maximum fold angle to F degrees, default is one degree"
    gdouble fold = M_PI/180.;
    //  fold = atof (optarg)*PI/180.;
    GtsStopFunc stop_func;
    gpointer stop_data = nullptr;
    if (edgeCountMax > 0) {
        stop_func = (GtsStopFunc) gts_coarsen_stop_number;
        stop_data = &edgeCountMax;
    } else {
        stop_func = (GtsStopFunc) gts_coarsen_stop_cost;
        stop_data = &costMax;
    }
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
             stop_func, stop_data, fold);

    auto newMesh(std::make_shared<geometry::Mesh>());
    fromGts (*newMesh, *s);
    gts_object_destroy (GTS_OBJECT (s));
    return newMesh;
}

gint getSubSurface (void *ee, gpointer *data)
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

    return {};
}

gint getSubSurfaces (GtsFace *f, gpointer *data)
{
    std::vector <std::vector <GtsSurface *>> *subSurfaces = (std::vector <std::vector <GtsSurface *>> *) data[0];
    std::function<math::Point2ll(double x, double y)> *getGridCell = (std::function<math::Point2ll(double, double)>* ) data[1];


    GtsVertex *v1, *v2, *v3;
    gts_triangle_vertices (GTS_TRIANGLE (f), &v1, &v2, &v3);

    gdouble midx = GTS_POINT (v1) -> x + GTS_POINT (v2)->x + GTS_POINT (v3) ->x;
    midx /= 3;
    gdouble midy = GTS_POINT (v1) -> y + GTS_POINT (v2)->y + GTS_POINT (v3) ->y;
    midy /= 3;

    auto targetCell = (*getGridCell)(midx, midy);
    GtsSurface *subSurface =(*subSurfaces)[targetCell(1)][targetCell(0)];
    gts_surface_add_face (subSurface, f);

    return {};
}

geometry::Mesh::pointer simplify_gts_in_grid(const geometry::Mesh &mesh
    , std::vector<std::vector <geometry::GridCell>> &gridCells
    , bool inParallel
    , std::function<math::Point2ll(double x, double y)> getGridCell)
{
    // GTS must be initialized at this point
    checkGtsInitialized();
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

    auto rows = gridCells.size();
    auto columns = gridCells[0].size(); // All rows have the same number of columns


    auto oneCellProcess ([&](long long  row, long long  column
                             , GtsSurface *subSurface) -> void
    {
        auto edgeMaxCount = gridCells[row][column].maxFaceCount * 3 / 2;
        // -1, because we are counting from 0!
        LOG (info3) << "Simplifying cell [" << row << "," << column << "]" << "/[" << rows - 1 << "," << columns - 1 << "] "
            << gridCells[row][column].extent
            << " to " << gridCells[row][column].maxFaceCount << " faces";
        gts_surface_coarsen (subSurface,
                cost_func, cost_data,
                coarsen_func, coarsen_data,
                stop_func, &edgeMaxCount, fold);

    });


    if (!inParallel) {
        for (unsigned int row = 0 ; row < rows; row+=1)
        for (unsigned int column = 0 ; column < columns; column+=1)
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
            gts_object_destroy (GTS_OBJECT (subSurface));
        }
    } else {
        LOG (info3) << "Creating empty subSurfaces";
        // Create a list of surfaces
        std::vector<std::vector <GtsSurface *>> subSurfaces(rows);
        for (unsigned int row = 0; row < rows ; ++row)
        for (unsigned int column = 0; column < columns ; ++column) {
            GtsSurface *subSurface = gts_surface_new (gts_surface_class (),
                    gts_face_class (),
                    gts_edge_class (),
                    gts_vertex_class ());
            subSurfaces[row].push_back (subSurface);
        }
        // Generate surfaces
        gpointer subParams[2];
        subParams[0] = &subSurfaces;
        subParams[1] = &getGridCell;
        LOG (info3) << "Filling subSurfaces";
        gts_surface_foreach_face (s, (GtsFunc) getSubSurfaces, &subParams);

        auto computeBunchInParallel = [&](long long cellRowStart
                                         , long long cellColumnStart
                                         , long long cellRowShift
                                         , long long cellColumnShift) -> void
        {
            LOG (info3)  << "Next bunch";
            // We cannot change a gts_surface (simplifying one cell of the mesh)
            // and at the same time iterate over gts_surface_verteces/faces/edges
            // (creating a subsurface for another cell)
            // The error:
            //     Gts:ERROR:surface.c:184:gts_surface_remove_face: assertion failed: (s->keep_faces == FALSE)
            // Solution: prepare all sub surfaces for this bunch in advance
#ifdef _MSC_VER
            UTILITY_OMP(parallel for schedule(dynamic, 1))
#else
            UTILITY_OMP(parallel for schedule(dynamic, 1) collapse(2))
#endif
            for (int64_t row = cellRowStart ; row < (int64_t)rows; row+=cellRowShift)
            for (unsigned int column = cellColumnStart ; column < columns;
                 column+=cellColumnShift)
            {
                GtsSurface *subSurface = subSurfaces[row][column];
                oneCellProcess (row, column, subSurface);
            }
            // Clean up to free some space
            for (unsigned int row = cellRowStart ; row < rows; row+=cellRowShift)
            for (unsigned int column = cellColumnStart ; column < columns;
                 column+=cellColumnShift) {
                GtsSurface *subSurface = subSurfaces[row][column];
                gts_object_destroy (GTS_OBJECT (subSurface));
            }
        };
        // These are guaranded to do not have intersection in 2-neghbourhoods
        // Because between every cells that are computed in parallel
        // there is at least one unsimplifyed cell, which has a loooot of triangles
        computeBunchInParallel (0, 0, 2, 2); // 1
        computeBunchInParallel (0, 1, 2, 4); // 2
        computeBunchInParallel (1, 0, 4, 4); // 3
        computeBunchInParallel (1, 1, 4, 4); // 4
        computeBunchInParallel (1, 2, 4, 4); // 5

        computeBunchInParallel (3, 0, 8, 4); // 6
        computeBunchInParallel (3, 1, 8, 4); // 7
        computeBunchInParallel (3, 2, 8, 4); // 8

        for (unsigned int row = 7; row < rows; row += 8) {
            computeBunchInParallel (row, 0, rows + 1, 4); // only this row  //6', 6'', 6''', ...
            computeBunchInParallel (row, 1, rows + 1, 4); // only this row  //7', 7'', 7''', ...
            computeBunchInParallel (row, 2, rows + 1, 4); // only this row  //8', 8'', 8''', ...
        }

        for (unsigned int row = 0; row < rows; ++row) {
            computeBunchInParallel (row, 3, rows + 1, 8); // only this row //9, 9', 9'', 9''', ...
        }

        // run the rest cells in sequence
        for (unsigned int row = 0 ; row < rows; row+=1)
        for (unsigned int column = 7 ; column < columns; column+=8)
        {
            GtsSurface *subSurface = subSurfaces[row][column];
            oneCellProcess (row, column, subSurface);
            gts_object_destroy (GTS_OBJECT (subSurface));
        }
    }

    auto newMesh(std::make_shared<geometry::Mesh>());
    fromGts (*newMesh, *s);
    gts_object_destroy (GTS_OBJECT (s));
    return newMesh;
}

} // namespace end

#include "gtsutils.hpp"
#include <dbglog/dbglog.hpp>

namespace geometry {

static void initializeGts()
{
    // Class system in GTS is not thread safe at compile time
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

static bool gtsInitialized = []{
    initializeGts();
    return true;
}();

void checkGtsInitialized() {
    if (!gtsInitialized) {
        LOGTHROW(err4, std::runtime_error) << "GTS has not been initialized";
    }
}

} // namespace geometry

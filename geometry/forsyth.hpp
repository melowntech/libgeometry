/**
 * Linear-Speed Vertex Cache Optimisation
 * https://tomforsyth1000.github.io/papers/fast_vert_cache_opt.html
 *
 * Used for optimal mesh serialization.
 *
 * \file forsyth.hpp
 * \author Vaclav Blazek <vaclav.blazek@melown.com>
 */

#ifndef geometry_forsyth_hpp_included
#define geometry_forsyth_hpp_included

#include <cstdint>

namespace geometry {

typedef std::uint16_t ForsythVertexIndexType;

bool forsythReorder(int *triOrder, const ForsythVertexIndexType *indices
                    , int nTriangles, int nVertices);

} // namespace geometry

#endif // geometry_forsyth_hpp_included

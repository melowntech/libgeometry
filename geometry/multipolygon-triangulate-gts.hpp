/**
 * Copyright (c) 2022 Melown Technologies SE
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
 *  @file triangulate.hpp
 *  @author Tomas Novak <tomas.novak@melowntech.com>
 *
 *  Multipolygon triangulation using CDT in GTS.
 */

#ifndef GEOMETRY_MULTIPOLYGON_TRIANGULATE_GTS_HPP_INCLUDED_
#define GEOMETRY_MULTIPOLYGON_TRIANGULATE_GTS_HPP_INCLUDED_

#include "math/geometry.hpp"
#include "utility/gccversion.hpp"

namespace geometry {

using PolygonIt = math::MultiPolygon::const_iterator;
using PointIt = math::Polygon::const_iterator;
using ItPair = std::pair<PolygonIt, PointIt>;

using TriangleItPair = std::array<ItPair, 3>;

/**
 * Triangulate a general multipolygon, return iterators pointing to input
 *
 * The function uses Constrained Delaunay triangulation implemented in GTS
 * library.
 *
 * @param[in] mpolygon input multipolygon (vector of polygons - first defines
 *                     outer boundary (CCW orientation), followed by inner
 *                     boundaries (CW orientation))
 * @return triangles with vertices defined by pairs of iterators - first
 *         pointing to a polygon in `mpolygon`, second to a vertex in polygon
 */
std::vector<TriangleItPair> multipolygonTriangulateGtsToIters(
    const math::MultiPolygon& mpolygon)
#ifndef GEOMETRY_HAS_GTS
    UTILITY_FUNCTION_ERROR(
        "Constrained Delaunay triangulation is available only when compiled "
        "with GTS.")
#endif
        ;

/**
 * Triangulate a general multipolygon
 *
 * The function uses Constrained Delaunay triangulation implemented in GTS
 * library.
 *
 * @param[in] mpolygon input multipolygon (vector of polygons - first defines
 *                     outer boundary (CCW orientation), followed by inner
 *                     boundaries (CW orientation))
 * @return triangles composing the polygon (without holes)
 */
math::Triangles2d multipolygonTriangulateGts(const math::MultiPolygon& mpolygon)
#ifndef GEOMETRY_HAS_GTS
    UTILITY_FUNCTION_ERROR(
        "Constrained Delaunay triangulation is available only when compiled "
        "with GTS.")
#endif
        ;

} // namespace geometry

#endif /* GEOMETRY_MULTIPOLYGON_TRIANGULATE_GTS_HPP_INCLUDED_ */

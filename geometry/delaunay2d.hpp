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
/**
 * @file delaunay2d.hpp
 * @author Jakub Cerveny <jakub.cerveny@melown.com>
 *
 * Computation of 2D Delaunay triangulation using the CGAL library.
 *
 * Note: you need to have CGAL configured (CGAL_FOUND in CMakeLists.txt)
 *       to use this module.
 */

#ifndef geometry_delaunay_hpp_included_
#define geometry_delaunay_hpp_included_

#include <array>
#include <vector>

#include "math/geometry_core.hpp"

#include "utility/gccversion.hpp"

namespace geometry {

typedef std::array<unsigned, 3> DTriangle;
typedef std::pair<unsigned, unsigned> DEdge;

/** Calculates the 2D Delaunay triangulation of a set of points.
 *  Returns a list of (finite) triangles. Each triangle indexes three
 *  points of the original set.
 */
std::vector<DTriangle>
    delaunayTriangulation2d(const math::Points2 &points)
#ifndef GEOMETRY_HAS_CGAL
    UTILITY_FUNCTION_ERROR("Delaunay triangulation is only available when "
                           "compiling with CGAL.")
#endif
    ;

/** Calculates the 2D constrained Delaunay triangulation of a set of points,
 *  where some of the edges are prescribed and will appear in the output.
 *  Returns a list of (finite) triangles. Each triangle references three
 *  points of the original set.
 */
/*std::vector<DTriangle>
    constrainedDelaunayTriangulation2d(
            const math::Points2 &points,
            const std::vector<DEdge> &constrained_edges)
#ifndef GEOMETRY_HAS_CGAL
    UTILITY_FUNCTION_ERROR("Constrained Delaunay triangulation is only "
                           "available when compiling with CGAL.")
#endif
    ;*/

} // namespace geometry

#endif // geometry_delaunay_hpp_included_

/**
 * @file delaunay2d.hpp
 * @author Jakub Cerveny <jakub.cerveny@ext.citationtech.net>
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

//! Calculates the 2D Delaunay triangulation of a set of points.
//! Returns a list of (finite) triangles. Each triangle references three
//! points of the original set.
//!
std::vector<DTriangle> delaunayTriangulation2d(const math::Points2 points)
#ifndef GEOMETRY_HAS_CGAL
    UTILITY_FUNCTION_ERROR("Delaunay triangulation is only available when"
                           " compiling with CGAL.")
#endif
    ;

} // namespace geometry

#endif // geometry_delaunay_hpp_included_

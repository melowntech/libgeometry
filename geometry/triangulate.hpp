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
 *  @file triangulate.hpp
 *  @author Jakub Cerveny <jakub.cerveny@melown.com>
 *
 *  Triangulation polygons.
 */

#ifndef geometry_triangulate_hpp_included_
#define geometry_triangulate_hpp_included_

#include "math/geometry.hpp"

namespace geometry {

/** Perform triangulation of a simple polygon by the ear clipping algorithm (see
 *  https://en.wikipedia.org/wiki/Polygon_triangulation#Ear_clipping_method)
 *  The polygon must be CCW oriented, must not self-intersect and must not
 *  contain holes. It may be nonconvex.
 */
math::Triangles2d simplePolyTriangulate(const math::Polygon &polygon);

/** Triangulate a general multipolygon (possibly with holes). The function first
 *  performs Delaunay triangulation of all points in the multipolygon and then
 *  returns triangles that lie inside.
 */
math::Triangles2d generalPolyTriangulate(const math::MultiPolygon &mpolygon);


} // namespace geometry

#endif // geometry_triangulate_hpp_included_

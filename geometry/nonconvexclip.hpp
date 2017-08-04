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
 *  @file nonconvexclip.hpp
 *  @author Jakub Cerveny <jakub.cerveny@melown.com>
 *
 *  Triangle clipping by a general nonconvex region.
 */

#ifndef geometry_nonconvexclip_hpp_included_
#define geometry_nonconvexclip_hpp_included_

#include "math/geometry.hpp"

#include <tuple>

namespace geometry {

typedef std::vector<math::Points2d> Region; // TODO: include from somewhere else

/** Clips a 3D triangle by a region in the XY plane. The result may be zero or
 *  more triangles covering the result of the boolean operation. The clip region
 *  may consist of multiple nonconvex (but simple) polygons. TODO: holes?
 */
math::Triangles3d clipTriangleNonconvex(const math::Triangle3d &tri,
                                        const Region &clipRegion);

/** Similar to clipTriangleNonconvex, but handles texture coordinates as well.
 */
std::tuple<math::Triangles3d, math::Triangles2d>
    clipTexturedTriangleNonconvex(const math::Triangle3d &tri,
                                  const math::Triangle2d &uv,
                                  const Region &clipRegion);

} // namespace geometry

#endif // geometry_nonconvexclip_hpp_included_

/**
 * Copyright (c) 2023 Melown Technologies SE
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
 *  @file polyop.hpp
 *  @author Jakub Cerveny <jakub.cerveny@melowntech.com>
 *
 *  Polygon operations, based on Boost.Geometry.
 *
 */

#ifndef geometry_polyop_hpp_included_
#define geometry_polyop_hpp_included_

#include "math/geometry.hpp"

namespace geometry {

/**
 * Calculate and return the intersection of two (multi-)polygons.
 */
math::MultiPolygon intersectPolygons(const math::MultiPolygon &mp1,
                                     const math::MultiPolygon &mp2);

/**
 * Calculate the difference "mp1 - mp2", i.e., return the shape containing
 * points from mp1 that are not contained in mp2.
 */
math::MultiPolygon subtractPolygons(const math::MultiPolygon &mp1,
                                    const math::MultiPolygon &mp2);

/**
 * Calculate and return the union of (multi-)polygons mp1 and mp2.
 */
math::MultiPolygon unitePolygons(const math::MultiPolygon &mp1,
                                 const math::MultiPolygon &mp2);

/**
 * Calculate the offset polygon (or buffer) of 'mpoly'.
 *
 * @param distance -- distance of the offset shape, can be both positive
 *   (inflate the polygon) or negative (reduce the polygon).
 * @param miterLimit -- limits the size of sharp "spikes".
 * @param epsSimplify -- simplification tolerance (0 = disable)
 */
math::MultiPolygon offsetPolygon(const math::MultiPolygon &mpoly,
                                 double distance, double miterLimit = 5,
                                 double epsSimplify = 1e-5);


} // namespace geometry

#endif // geometry_polyop_hpp_included_

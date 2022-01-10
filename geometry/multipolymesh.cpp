/**
 * Copyright (c) 2021 Melown Technologies SE
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

#include "multipolymesh.hpp"

#include "math/math.hpp"

#include <boost/numeric/ublas/matrix.hpp>

namespace geometry
{
namespace ublas = boost::numeric::ublas;

FacePlaneCrs::FacePlaneCrs(const math::Point3& p1,
                           const math::Point3& p2,
                           const math::Point3& p3)
{
    // get base vectors
    math::Point3 n1 = math::normalize(p2 - p1);
    normal_ = math::normalize(math::crossProduct(p2 - p1, p3 - p1));
    math::Point3 n2 = math::normalize(math::crossProduct(normal_, n1));

    p2g_ = math::identity4();

    auto col1 = ublas::column(p2g_, 0);
    auto col2 = ublas::column(p2g_, 1);
    auto col3 = ublas::column(p2g_, 2);
    auto col4 = ublas::column(p2g_, 3);
    ublas::subrange(col1, 0, 3) = n1;
    ublas::subrange(col2, 0, 3) = n2;
    ublas::subrange(col3, 0, 3) = normal_;
    ublas::subrange(col4, 0, 3) = p1;

    g2p_ = math::matrixInvert(p2g_);
}

} // namespace geometry

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
 *  @file polyclip.cpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Polygon clipping stuff
 *
 *  2013-01-03 (vasek)    created
 */

#include "dbglog/dbglog.hpp"

#include "polygon.hpp"


namespace geometry { namespace detail {

namespace cohen_hodgman {

enum { Right, Left, Upper, Lower };

class Clipper {
public:
    Clipper(const math::Extents2 &ext) : ext(ext) {}

    math::Points2 clip(const math::Points2 &polygon) {
        return clipHalfspace<Lower>
            (clipHalfspace<Left>
             (clipHalfspace<Upper>
              (clipHalfspace<Right>(polygon))));
    }

private:
    const math::Extents2 ext;

    template <int Halfspace>
    bool inside(const math::Point2 &p);

    template <int Halfspace>
    math::Point2 intersect(const math::Point2 &p0, const math::Point2 &p1);

    template <int Halfspace>
    math::Points2 clipHalfspace(const math::Points2 &polygon);
};

template <>
bool Clipper::inside<Right>(const math::Point2 &p)
{
    return p(0) >= ext.ll(0);
}

template <>
bool Clipper::inside<Left>(const math::Point2 &p)
{
    return p(0) <= ext.ur(0);
}

template <>
bool Clipper::inside<Upper>(const math::Point2 &p)
{
    return p(1) >= ext.ll(1);
}

template <>
bool Clipper::inside<Lower>(const math::Point2 &p)
{
    return p(1) <= ext.ur(1);
}

template <>
math::Point2 Clipper::intersect<Right>(const math::Point2 &p0
                                       , const math::Point2 &p1)
{
    return {ext.ll(0)
            , p0(1) + (p1(1) - p0(1)) * (ext.ll(0) - p0(0))
            / (p1(0) - p0(0))};
}

template <>
math::Point2 Clipper::intersect<Left>(const math::Point2 &p0
                                      , const math::Point2 &p1)
{
    return {ext.ur(0)
            , p0(1) + (p1(1) - p0(1)) * (ext.ur(0) - p0(0))
            / (p1(0) - p0(0))};
}

template <>
math::Point2 Clipper::intersect<Upper>(const math::Point2 &p0
                                       , const math::Point2 &p1)
{
    return { p0(0) + (p1(0) - p0(0)) * (ext.ll(1) - p0(1))
            / (p1(1) - p0(1))
            , ext.ll(1) };
}

template <>
math::Point2 Clipper::intersect<Lower>(const math::Point2 &p0
                                       , const math::Point2 &p1)
{
    return { p0(0) + (p1(0) - p0(0)) * (ext.ur(1) - p0(1))
            / (p1(1) - p0(1))
            , ext.ur(1) };
}

template <int Halfspace>
math::Points2 Clipper::clipHalfspace(const math::Points2 &polygon)
{
    if (polygon.empty()) {
        return {};
    }

    math::Points2 res;
    res.reserve(2 * polygon.size());

    // previous point
    auto prev(polygon.back());
    bool prevInside(inside<Halfspace>(prev));

    for (const auto &p : polygon) {
        bool pInside(inside<Halfspace>(p));
        // "logical xor"
        if (pInside != prevInside) {
            // line segment crosses halfspace boundary

            // add intersection
            res.push_back(intersect<Halfspace>(prev, p));
        }

        if (pInside) {
            // add p as well
            res.push_back(p);
        }

        // prepare for next round
        prev = p;
        prevInside = pInside;
    }

    return res;
}

} // namespace cohen_hodgman

math::Points2 clip(const math::Viewport2f &viewport
                   , const math::Points2 &polygon)
{
    return cohen_hodgman::Clipper(extents(viewport)).clip(polygon);
}

} } // namespace geometry::detail

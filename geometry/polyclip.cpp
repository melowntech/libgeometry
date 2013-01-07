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

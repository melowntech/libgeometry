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
 *  @file geos-geometry-convert.cpp
 *
 *  Conversions from and to GEOS geometry.
 *
 */

#include "geos-geometry-convert.hpp"

#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/CoordinateArraySequence.h>

#include <cassert>

namespace geometry {

using namespace geos::geom;
std::unique_ptr<Polygon> convert2geos(const math::Triangle2d &tri) 
{
    // Note: ignoring orientation
    auto g = GeometryFactory::create();
    
    std::unique_ptr<CoordinateArraySequence> seq(new CoordinateArraySequence(4,2));
    seq->setAt({tri[0][0], tri[0][1]}, 0);
    seq->setAt({tri[1][0], tri[1][1]}, 1);
    seq->setAt({tri[2][0], tri[2][1]}, 2);
    seq->setAt({tri[0][0], tri[0][1]}, 3);

    std::unique_ptr<Polygon> poly = g->createPolygon(g->createLinearRing(std::move(seq)));

    return poly;
}

math::Triangle2d convert2math(const geos::geom::Polygon *tri) 
{
    // Note: ignoring orientation
    assert(tri->getExteriorRing()->getNumPoints() == 4 && tri->getNumInteriorRing() == 0);
    math::Triangle2d t{
       math::Point2d{tri->getExteriorRing()->getCoordinateN(0).x, tri->getExteriorRing()->getCoordinateN(0).y},
       math::Point2d{tri->getExteriorRing()->getCoordinateN(1).x, tri->getExteriorRing()->getCoordinateN(1).y},
       math::Point2d{tri->getExteriorRing()->getCoordinateN(2).x, tri->getExteriorRing()->getCoordinateN(2).y}
    };
    return t;
}

// Note:
//
//   TODO: math::MultiPolygon should be defined better without ambiguities !!!
//   math::Polygon            - single CCW ring, not closed
//   math::MultiPolygon       - multiple rings, holes CW
//
//   geos::geom::Polygon      - single shell (closed, un-oriented), any number of holes (closed, un-oriented)
//   geos::geom::MultiPolygon - multiple polygons
//

std::unique_ptr<MultiPolygon> convert2geos(const math::MultiPolygon &mpoly) 
{
    // TODO: find proper solution for conversion to geos and
    // resolve ambiguities in math::MultiPolygon
    // --------------------------------------------------

    // Current implementation assumes that math::MultiPolygon
    // is a sequence of polygon rings (CCW, not closed ring) and 
    // holes (CW, not closed ring) and individual polygons are 
    // separated by the rings in the sequence.

    // Examples:
    // math::MultiPolygon
    //  {CW, CW, CW, CCW, CW, CW, CCW, CCW}
    // geos::geom::MultiPolygon
    //  {{1 ring, 3 holes}, {1 ring, 2 holes}, {1 ring}}
    //
    // math::MultiPolygon
    //  {CCW, CCW, CW, CW, CW, CW, CW, CCW, CCW}
    // geos::geom::MultiPolygon
    //  {{1 ring, 0 holes}, {1 ring, 5 holes}, {1 ring, 0 holes}, {1 ring, 0 holes}}

    auto g = GeometryFactory::create();
    std::vector<std::unique_ptr<Polygon>> polygons;
    std::vector<std::unique_ptr<LinearRing>> rings, holes;
    for (auto ring : mpoly) {
        if (ring.size()) {
            std::unique_ptr<CoordinateArraySequence> seq(new CoordinateArraySequence(ring.size() + 1, 2));
            seq->setAt({ring.at(0)[0], ring.at(0)[1]}, 0);
            double cw = 0.0;
            for (std::size_t i = 1; i <= ring.size(); i++) {
                auto const & prev = ring.at((i - 1) % ring.size());
                auto const & current = ring.at(i % ring.size());
                seq->setAt({current[0], current[1]}, i);
                cw += (current[0] - prev[0]) * (current[1] + prev[1]);
            }
            if (cw > 0.0) {
                holes.emplace_back(g->createLinearRing(std::move(seq)));
            } else {
                if (rings.size()) {
                    polygons.emplace_back(g->createPolygon(std::move(rings.back()), std::move(holes)));
                    rings.clear();
                    holes.clear();
                }
                rings.emplace_back(g->createLinearRing(std::move(seq)));
            }
        }
    }
    if (rings.size()) {
        polygons.emplace_back(g->createPolygon(std::move(rings.back()), std::move(holes)));
        rings.clear();
        holes.clear();
    }

    return g->createMultiPolygon(std::move(polygons));
}


math::MultiPolygon convert2math(const geos::geom::Geometry *g)
{
    math::MultiPolygon result;

    auto addRing = [&](const geos::geom::LinearRing *ring, bool cw)
    {
        const geos::geom::CoordinateSequence* seq = ring->getCoordinatesRO();
        if (seq->getSize()) {
            math::Polygon r;
            // check ring orientation
            double ringCW = 0.0;
            for (std::size_t i = 1; i < seq->getSize(); i++) {
                ringCW += (seq->getAt(i).x - seq->getAt(i - 1).x) * (seq->getAt(i).y + seq->getAt(i - 1).y);
            }
            if ((ringCW > 0.0) == cw) {
                // desired orientation matches
                // add current ring (without last point, we want open ring)
                for (std::size_t i = 0; i <= seq->getSize() - 2; i++) {
                    r.emplace_back(math::Point2d{seq->getAt(i).x, seq->getAt(i).y});
                }
            } else {
                // desired orientation does not match
                // add reverse of current ring (without first point, we want open ring)
                for (std::size_t i = seq->getSize() - 1; i >= 1; i--) {
                    r.emplace_back(math::Point2d{seq->getAt(i).x, seq->getAt(i).y});
                }
            }
            result.emplace_back(r);
        }
    };

    auto addPoly = [&](const geos::geom::Polygon* poly)
    {
        const geos::geom::LinearRing *shell = poly->getExteriorRing();
        addRing(shell, false);
        for (std::size_t i = 0; i < poly->getNumInteriorRing(); i++) {
            const geos::geom::LinearRing *hole = poly->getInteriorRingN(i);
            addRing(hole, true);
        }
    };

    if (g->getGeometryTypeId() == geos::geom::GeometryTypeId::GEOS_POLYGON) {
        const geos::geom::Polygon* poly = static_cast<const geos::geom::Polygon*>(g);
        addPoly(poly);
    }
    else if (g->getGeometryTypeId() == geos::geom::GeometryTypeId::GEOS_MULTIPOLYGON) {
        const geos::geom::MultiPolygon* multiPoly = static_cast<const geos::geom::MultiPolygon*>(g);
        for (std::size_t i = 0; i < multiPoly->getNumGeometries(); i++) {
            const geos::geom::Polygon* poly = multiPoly->getGeometryN(i);
            addPoly(poly);
        }
    } else {
        // TODO: warning?
    }

    return result;
}

} // namespace geometry

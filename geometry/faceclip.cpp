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
#include <cmath>

#include "faceclip.hpp"

namespace geometry { namespace opencv {

typedef ClipTriangle::Point Point;

namespace {

inline double signedDistance(const Point &point, const ClipPlane &plane)
{
    return point.dot(plane.normal) + plane.d;
}


Point intersection(const Point &p1, const Point &p2,
                   const ClipPlane &plane, double& t)
{
    double dot1 = p1.dot(plane.normal);
    double dot2 = p2.dot(plane.normal);
    double den = dot1 - dot2;

    // line parallel with plane, return the midpoint
    if (std::abs(den) < 1e-10) {
        t = 0.5;
        return (p1 + p2) * t;
    }

    t = (dot1 + plane.d) / den;
    return (1.0 - t)*p1 + t*p2;
}

} // namespace


ClipTriangle::list clipTriangles(const ClipTriangle::list &triangles,
                                 const ClipPlane &plane)
{
    ClipTriangle::list result;

    for (const auto &tri : triangles)
    {
        bool positive[3] = {
            signedDistance(tri.pos[0], plane) >= 0,
            signedDistance(tri.pos[1], plane) >= 0,
            signedDistance(tri.pos[2], plane) >= 0
        };

        int count = 0;
        for (int i = 0; i < 3; i++) {
            if (positive[i]) count++;
        }

        // triangle completely on negative side - do nothing
        if (count == 0) continue;

        // trinagle completely on positive side - copy to result
        if (count == 3) {
            result.push_back(tri);
            continue;
        }

        int a, b, c;
        double t;

        // case 1: one vertex on positive side, just adjust the other two
        if (count == 1)
        {
            if (positive[0]) a = 0, b = 1, c = 2;
            else if (positive[1]) a = 1, b = 2, c = 0;
            else a = 2, b = 0, c = 1;

            auto x1pos(intersection(tri.pos[a], tri.pos[b], plane, t));
            auto x1uv((1.0 - t)*tri.uv[a] + t*tri.uv[b]);

            auto x2pos(intersection(tri.pos[c], tri.pos[a], plane, t));
            auto x2uv((1.0 - t)*tri.uv[c] + t*tri.uv[a]);

            result.emplace_back(tri.id1, tri.id2,
                                tri.pos[a], x1pos, x2pos,
                                tri.uv[a],  x1uv,  x2uv);
        }
        // case 2: two vertices on positive side, adjust triangle and add one more
        else
        {
            if (!positive[0]) a = 0, b = 1, c = 2;
            else if (!positive[1]) a = 1, b = 2, c = 0;
            else a = 2, b = 0, c = 1;

            auto x1pos(intersection(tri.pos[a], tri.pos[b], plane, t));
            auto x1uv((1.0 - t)*tri.uv[a] + t*tri.uv[b]);

            auto x2pos(intersection(tri.pos[c], tri.pos[a], plane, t));
            auto x2uv((1.0 - t)*tri.uv[c] + t*tri.uv[a]);

            result.emplace_back(tri.id1, tri.id2,
                                x1pos, tri.pos[b], tri.pos[c],
                                x1uv,  tri.uv[b],  tri.uv[c]);

            result.emplace_back(tri.id1, tri.id2,
                                x1pos, tri.pos[c], x2pos,
                                x1uv,  tri.uv[c],  x2uv);
        }
    }

    return result;
}

} } // namespace geometry::opencv

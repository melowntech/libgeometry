/**
 * Copyright (c) 2022 Melown Technologies SE
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

#include "nonconvexclip.hpp"

#include "meshop.hpp"

namespace geometry {

namespace {

// 1 micropixel
const double CompareEpsilon(1e-6);

inline bool doubleLess(double a, double b) {
    if (std::abs(a - b) <= CompareEpsilon) { return false; }
    return a < b;
}

struct PointCompare {
    bool operator()(const math::Point3d &a, const math::Point3d &b) const {
        if (doubleLess(a(0), b(0))) { return true; }
        if (doubleLess(b(0), a(0))) { return false; }

        if (doubleLess(a(1), b(1))) { return true; }
        if (doubleLess(b(1), a(1))) { return false; }

        return doubleLess(a(2), b(2));
    }

    bool operator()(const math::Point2d &a, const math::Point2d &b) const {
        if (doubleLess(a(0), b(0))) { return true; }
        if (doubleLess(b(0), a(0))) { return false; }

        return doubleLess(a(1), b(1));
    }
};

inline int which(const math::Triangle3d &t, const math::Point3 &p)
{
    if (p == t[0]) { return 0; }
    if (p == t[1]) { return 1; }
    if (p == t[2]) { return 2; }
    return -1;
}

inline int which(const math::Triangle2d &t, const math::Point2 &p)
{
    if (p == t[0]) { return 0; }
    if (p == t[1]) { return 1; }
    return -1;
}

} // namespace

Mesh clip(const Mesh &mesh, const math::MultiPolygon &clipRegion)
{
    geometry::Mesh om;
    using VertexMap = std::vector<int>;
    VertexMap vertexMap(mesh.vertices.size(), -1);
    VertexMap uvMap(mesh.tCoords.size(), -1);

    const auto addVertexIndex([&](std::size_t v) -> std::size_t
    {
        auto &mapping(vertexMap[v]);
        if (mapping < 0) {
            mapping = om.vertices.size();
            om.vertices.push_back(mesh.vertices[v]);
        }
        return mapping;
    });

    const auto addUvIndex([&](std::size_t v) -> std::size_t
    {
        auto &mapping(uvMap[v]);
        if (mapping < 0) {
            mapping = om.tCoords.size();
            om.tCoords.push_back(mesh.tCoords[v]);
        }
        return mapping;
    });

    // NB: here, we expect that original vertices are kept as-is therefore we
    // can map them to original vertices
    using SyntheticVertices = std::map<math::Point3d, int, PointCompare>;
    SyntheticVertices syntheticVertices;

    const auto addSyntheticVertex([&](const math::Point3d &p) -> std::size_t
    {
        auto fsyntheticVertices(syntheticVertices.find(p));
        if (fsyntheticVertices != syntheticVertices.end()) {
            return fsyntheticVertices->second;
        }

        const auto index(om.vertices.size());
        syntheticVertices.insert(
            SyntheticVertices::value_type(p, static_cast<int>(index)));
        om.vertices.push_back(p);
        return index;
    });

    using SyntheticUv = std::map<math::Point2d, int, PointCompare>;
    SyntheticUv syntheticUv;

    const auto addSyntheticUv([&](const math::Point2d &p) -> std::size_t
    {
        auto fsyntheticUv(syntheticUv.find(p));
        if (fsyntheticUv != syntheticUv.end()) {
            return fsyntheticUv->second;
        }

        const auto index(om.tCoords.size());
        syntheticUv.insert(SyntheticUv::value_type(p, index));
        om.tCoords.push_back(p);
        return index;
    });

    const bool textured(!mesh.tCoords.empty());

    for (const auto &face : mesh.faces) {
        // clip face

        const std::array<std::size_t, 3>
            fIndices{{ face.a, face.b, face.c }};

        const math::Triangle3d f{{ mesh.vertices[face.a]
                                   , mesh.vertices[face.b]
                                   , mesh.vertices[face.c] }};

        const auto addVertex([&](const math::Point3d &p) -> std::size_t
        {
            auto index(which(f, p));
            if (index < 0) { return addSyntheticVertex(p); }
            return addVertexIndex(fIndices[index]);
        });

        if (!textured) {
            // non-textured mesh
            const auto t3d(geometry::clipTriangleNonconvex(f, clipRegion));

            for (const auto &mt : t3d) {
                om.faces.emplace_back
                    (static_cast<geometry::Face::index_type>(addVertex(mt[0]))
                    , static_cast<geometry::Face::index_type>(addVertex(mt[1]))
                    , static_cast<geometry::Face::index_type>(
                        addVertex(mt[2])));
            }
            continue;
        }

        // textured mesh

        const std::array<std::size_t, 3>
            uvIndices{{ face.ta, face.tb, face.tc }};

        const math::Triangle2d t{{ mesh.tCoords[face.ta]
                                   , mesh.tCoords[face.tb]
                                   , mesh.tCoords[face.tc] }};

        const auto addUv([&](const math::Point2d &p) -> std::size_t
        {
            auto index(which(t, p));
            if (index < 0) { return addSyntheticUv(p); }
            return addUvIndex(uvIndices[index]);
        });

        math::Triangles3d t3d;
        math::Triangles2d t2d;
        std::tie(t3d, t2d)
            = geometry::clipTexturedTriangleNonconvex(f, t, clipRegion);

        auto it2d(t2d.begin());
        for (const auto &mt : t3d) {
            const auto &tt(*it2d++);
            om.faces.emplace_back
                (addVertex(mt[0]), addVertex(mt[1]), addVertex(mt[2])
                 , addUv(tt[0]), addUv(tt[1]), addUv(tt[2])
                 , face.imageId);
        }
    }

    return om;
}

} // namespace geometry

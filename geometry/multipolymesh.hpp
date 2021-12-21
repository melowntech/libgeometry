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
/**
 * @file multipolymesh.hpp
 * @author Tomas Novak <tomas.novak@melowntech.com>
 *
 * 3D mesh with multipolygonal faces (with holes)
 *
 */

#ifndef geometry_multipolymesh_hpp_included_
#define geometry_multipolymesh_hpp_included_

#include "geometry/mesh.hpp"
#include "geometry/multipolygon-triangulate-gts.hpp"
#include "geometry/polygon.hpp"
#include "math/geometry_core.hpp"
#include "math/transform.hpp"

namespace geometry
{
/**
 * Class specifying a 2D coordinate system of a planar face
 */
class FacePlaneCrs
{
    math::Matrix4 p2g_;
    math::Matrix4 g2p_;
    math::Point3 normal_;

public:
    /**
     * Initalizes the transformations to/from 2D crs given 3 points in plane.
     *
     * The points must not be collinear and must go in CCW direction along the
     * face boundary.
     */
    FacePlaneCrs(const math::Point3& p1,
                 const math::Point3& p2,
                 const math::Point3& p3);


    /**
     * Transforms a 3D point from global crs to the face 2D crs
     */
    inline math::Point2 to2D(const math::Point3& pt) const
    {
        math::Point3 pt2 = math::transform(g2p_, pt);
        return math::Point2(pt2(0), pt2(1));
    }

    /**
     * Transforms a 2D point from face 2D crs to the global crs
     */
    inline math::Point3 to3D(const math::Point2& pt) const
    {
        math::Point3 pt3(pt(0), pt(1), 0.0);
        return math::transform(p2g_, pt3);
    }

    /**
     * Returns the normal vector of plane
     */
    inline const math::Point3& normal() const { return normal_; }
};

/**
 * Face of multipolygonal mesh
 *
 * First chain - outer boundary (CCW) - is followed by holes (CW).
 */
using MultiPolyFace = std::vector<std::vector<std::size_t>>;

/**
 * Mesh with multipolygonal faces (i.e. polygons with holes)
 *
 * Optionally stores semantic labels and coordiante systems of individual faces.
 */
template <typename FaceLabels>
struct MultiPolyMesh
{
    /**
     * Vertices of mesh
     */
    math::Points3 vertices;

    /**
     * Multipolygonal faces of mesh
     */
    std::vector<MultiPolyFace> faces;

    /**
     * Optional labels of faces (may be empty)
     */
    std::vector<FaceLabels> faceLabels;

    /**
     * Optional 2D coordinate systems of faces (may be empty), also specify
     * normals of faces
     */
    std::vector<FacePlaneCrs> faceCrs;

    /**
     * Creates coordinate system object for a face
     *
     * Selects three suitable vertices (non-collinear) and constructs a
     * FacePlaneCrs
     *
     * @param[in] faceIdx face index
     * @returns the coordinate system object
     */
    FacePlaneCrs getFaceCrs(std::size_t faceIdx) const;

    /**
     * Creates coordinate systems for all faces and stores them to `faceCrs`.
     *
     * `faceCrs` is cleared in the process
     */
    void initalizeFaceCrs();

    /**
     * Triangulate polygonal faces
     *
     * If `faceCrs` is empty, calls getFaceCrs()
     *
     * @param[out] triFaceLabels optional labels of faces in triangular mesh
     *                           (generated if `faceLabels` is not empty)
     * @returns triangular mesh
     */
    Mesh triangulateFaces(std::vector<FaceLabels>* const triFaceLabels
                          = nullptr) const;
};

// Implementation following

namespace detail
{
constexpr double eps = 1e-9;
}

template <typename FaceLabels>
FacePlaneCrs MultiPolyMesh<FaceLabels>::getFaceCrs(std::size_t faceIdx) const
{
    auto& ob = faces[faceIdx][0];
    math::Point3 p1 = vertices[ob[0]];

    std::size_t i = 1;
    math::Point3 p2 = vertices[ob[i]];

    while (math::length(p2 - p1) < detail::eps)
    {
        ++i;
        if (i >= ob.size())
        {
            LOGTHROW(err4, std::runtime_error)
                << "Cannot find another well-distant point";
        }
        p2 = vertices[ob[i]];
    }

    math::Line3 l(p1, p2 - p1);

    math::Point3 p3 = vertices[ob[i]];

    while (math::pointLineDistance(p3, l) < detail::eps)
    {
        ++i;

        if (i >= ob.size())
        {
            LOGTHROW(err4, std::runtime_error)
                << "Cannot find another non-collinear point";
        }

        p3 = vertices[ob[i]];
    }

    FacePlaneCrs crs(p1, p2, p3);

    // check that the outer boundary has positive area when transformed
    math::Points2 pts2D;
    for (auto& v : ob)
    {
        pts2D.push_back(crs.to2D(vertices[v]));
    }
    if (area(pts2D) > 0) { return crs; }
    return FacePlaneCrs(p1, p3, p2);
}

template <typename FaceLabels>
void MultiPolyMesh<FaceLabels>::initalizeFaceCrs()
{
    faceCrs.clear();
    faceCrs.reserve(faces.size());

    for (std::size_t fIdx = 0; fIdx < faces.size(); ++fIdx)
    {
        faceCrs.push_back(getFaceCrs(fIdx));
    }
}

/// Helper - convert iterator pair in multipolygon to index of vertex in
/// MultiPolyMesh
inline std::size_t itPair2Vertex(const ItPair& itPair,
                                 const math::MultiPolygon& poly2D,
                                 const MultiPolyFace& polyFace)
{
    std::size_t polyIdx = std::distance(poly2D.begin(), itPair.first);
    std::size_t ptIdx = std::distance(poly2D[polyIdx].begin(), itPair.second);
    return polyFace[polyIdx][ptIdx];
}

template <typename FaceLabels>
Mesh MultiPolyMesh<FaceLabels>::triangulateFaces(
    std::vector<FaceLabels>* const triFaceLabels) const
{
    Mesh triMesh;
    if (triFaceLabels) { triFaceLabels->clear(); }
    triMesh.vertices = vertices; // copy vertices

    // triangulate polygons
    for (std::size_t fIdx = 0; fIdx < faces.size(); fIdx++)
    {
        auto& mpFace = faces[fIdx];

        FacePlaneCrs crsFace(faceCrs.empty() ? getFaceCrs(fIdx)
                                             : faceCrs[fIdx]);

        // get multipolygon in 2D
        math::MultiPolygon poly2D(mpFace.size());

        for (std::size_t i = 0; i < mpFace.size(); ++i)
        {
            poly2D[i].reserve(mpFace[i].size());
            std::transform(mpFace[i].begin(),
                           mpFace[i].end(),
                           std::back_inserter(poly2D[i]),
                           [this, &crsFace](const std::size_t vIdx) {
                               return crsFace.to2D(vertices[vIdx]);
                           });
        }

        std::vector<TriangleItPair> triangles
            = multipolygonTriangulateGts(poly2D);

        for (auto& t : triangles)
        {
            std::size_t v1 = itPair2Vertex(t[0], poly2D, mpFace);
            std::size_t v2 = itPair2Vertex(t[1], poly2D, mpFace);
            std::size_t v3 = itPair2Vertex(t[2], poly2D, mpFace);
            triMesh.addFace(v1, v2, v3);
            if (faceLabels.size() && triFaceLabels)
            {
                triFaceLabels->push_back(faceLabels[fIdx]);
            }
        }
    }
    return triMesh;
}

} // namespace geometry

#endif /* geometry_multipolymesh_hpp_included_ */

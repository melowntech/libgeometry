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

#include "mesh-polygonization.hpp"

#include <queue>
#include <cmath>
#include <map>

#include "math/math.hpp"
#include "math/geometry.hpp"
#include "utility/small_set.hpp"

#include "polygon.hpp"

namespace geometry
{

namespace
{

/// Multipolygonal face (i.e. polygonal face with holes) with vertices in pmesh
using VhMpolyFace = std::vector<std::vector<OMPolyMesh::VertexHandle>>;

inline math::Point3 fromOM(const OMPolyMesh::Point& p)
{
    return math::Point3(p[0], p[1], p[2]);
}

inline bool isFinite(const OMPolyMesh::Normal& n)
{
    return std::isfinite(n[0]) && std::isfinite(n[1]) && std::isfinite(n[2]);
}


/// Returns true if the halfedge is on region boundary. Optionally checks for
/// specific left region.
inline bool isBoundaryHalfedge(const OMPolyMesh& pmesh,
                               const OMPolyMesh::HalfedgeHandle& heh,
                               const OMFacePropInt& faceRegionProp,
                               int checkLeftRegion = -1)
{
    auto leftR { pmesh.property(faceRegionProp, pmesh.face_handle(heh)) };
    if (checkLeftRegion >= 0)
    {
        if (leftR != checkLeftRegion) { return false; }
    }

    auto op { pmesh.opposite_halfedge_handle(heh) };
    if (pmesh.is_boundary(op))
    {
        LOG(warn2) << "Mesh is not watertight.";
        return true;
    }
    auto rightR { pmesh.property(faceRegionProp, pmesh.face_handle(op)) };
    return leftR != rightR;
}

/// Returns true if the vertex is adjacent to more than two regions
inline bool isImportantVertex(const OMPolyMesh& pmesh,
                              const OMPolyMesh::VertexHandle& vh,
                              const OMFacePropInt& faceRegionProp)
{
    int r1 { -1 };
    int r2 { -1 };

    // circulate around the current vertex
    for (auto vfit = pmesh.cvf_iter(vh); vfit.is_valid(); ++vfit)
    {
        int region = pmesh.property(faceRegionProp, *vfit);
        if (region == r1 || region == r2) { continue; }
        if (r1 == -1)
        {
            r1 = region;
            continue;
        }
        if (r2 == -1)
        {
            r2 = region;
            continue;
        }
        return true;
    }
    return false;
}

/// Find the next halfedge that points to a face of the same region
inline OMPolyMesh::HalfedgeHandle
    getNextRegionHalfedge(const OMPolyMesh& pmesh,
                          const OMPolyMesh::HalfedgeHandle& hehPrev,
                          const OMFacePropInt& faceRegionProp,
                          bool searchCCW = false)
{
    auto rIdx { pmesh.property(faceRegionProp, pmesh.face_handle(hehPrev)) };

    // iterate over outcoming edges
    auto hehStart { pmesh.next_halfedge_handle(hehPrev) };
    if (searchCCW)
    {
        hehStart = pmesh.opposite_halfedge_handle(pmesh.prev_halfedge_handle(
            pmesh.opposite_halfedge_handle(hehPrev)));
    }

    auto heh { hehStart };
    do
    {
        if (isBoundaryHalfedge(pmesh, heh, faceRegionProp, rIdx))
        {
            return heh;
        }

        if (!searchCCW)
        {
            heh = pmesh.next_halfedge_handle(
                pmesh.opposite_halfedge_handle(heh));
        }
        else
        {
            heh = pmesh.opposite_halfedge_handle(
                pmesh.prev_halfedge_handle(heh));
        }
    }
    while (heh != hehStart);

    LOGTHROW(err4, std::runtime_error)
        << "Cannot find the next halfedge on the region boundary (something is "
           "very wrong).";
    throw;
}

/// Traverse a region boundary, recording vertices and adding visited halfedges
/// to `visitedHalfedges`. Check for
void traverseRegionBoundary(
    const OMPolyMesh& pmesh,
    const OMPolyMesh::HalfedgeHandle& hehStart,
    const OMFacePropInt& faceRegionProp,
    std::vector<std::vector<OMPolyMesh::VertexHandle>>& boundaryVertices,
    utility::small_set<OMPolyMesh::HalfedgeHandle>& visitedHalfedges,
    bool knownIsInnerBoundary = false)
{
    utility::small_set<OMPolyMesh::VertexHandle> visitedVertices;
    auto heh { hehStart };

    auto duplicateVertex { false };
    std::vector<OMPolyMesh::VertexHandle> boundary; // add new boundary
    do
    {
        visitedHalfedges.insert(heh);
        auto v { pmesh.from_vertex_handle(heh) };

        if (visitedVertices.count(v))
        {
            // visiting same vertex twice - holes meeting at vertex
            duplicateVertex = true;
            break;
        }

        visitedVertices.insert(v);
        if (isImportantVertex(pmesh, v, faceRegionProp))
        {
            boundary.push_back(v);
        }
        heh = getNextRegionHalfedge(pmesh,
                                    heh,
                                    faceRegionProp,
                                    knownIsInnerBoundary);
    }
    while (heh != hehStart);

    if (duplicateVertex && knownIsInnerBoundary)
    {
        LOGTHROW(err4, std::runtime_error)
            << "Encountering duplicate vertex when traversing known inner "
               "boundary";
    }

    if (duplicateVertex)
    {
        // restart traversal from all outgoing edges of the duplicate vertex
        auto v { pmesh.from_vertex_handle(heh) };
        auto rIdx { pmesh.property(faceRegionProp,
                                   pmesh.face_handle(hehStart)) };

        for (auto it { pmesh.cvoh_begin(v) }; it.is_valid(); ++it)
        {
            if (!isBoundaryHalfedge(pmesh, *it, faceRegionProp, rIdx))
            {
                continue;
            }
            traverseRegionBoundary(pmesh,
                                   *it,
                                   faceRegionProp,
                                   boundaryVertices,
                                   visitedHalfedges,
                                   true);
        }
    }
    else
    {
        boundaryVertices.push_back(std::move(boundary));
    }
}

/// Returns area of polygon defined by vertices of `boundary` transformed to 2D
/// by `tf`.
double boundaryArea(const std::vector<OMPolyMesh::VertexHandle>& boundary,
                    const OMPolyMesh& pmesh,
                    const math::Matrix4& glob2plane)
{
    math::Points2 pts;
    pts.reserve(boundary.size());
    for (const auto& v : boundary)
    {
        math::Point3 pt3(math::transform(glob2plane, fromOM(pmesh.point(v))));
        pts.emplace_back(pt3(0), pt3(1));
    }
    return geometry::area(pts);
}

/// Traverse connected faces of one region and create a multipolygonal face
VhMpolyFace traverseConnectedFacesOfRegion(
    const OMPolyMesh::HalfedgeHandle& hehStart,
    const math::Matrix4& glob2plane,
    const OMPolyMesh& pmesh,
    const OMFacePropInt& faceRegionProp,
    utility::small_set<OMPolyMesh::HalfedgeHandle>& traversedHalfedges)
{
    // get main info
    int rIdx = pmesh.property(faceRegionProp, pmesh.face_handle(hehStart));

    // Traverse whole connected component of the region
    utility::small_set<OMPolyMesh::HalfedgeHandle> enqueuedHalfedges;
    std::queue<OMPolyMesh::HalfedgeHandle> halfedgeQ;
    halfedgeQ.push(hehStart);
    enqueuedHalfedges.insert(hehStart);

    int outerBoundaryNum { 0 };
    VhMpolyFace multipolygon;

    while (!halfedgeQ.empty())
    {
        auto heh { halfedgeQ.front() };
        halfedgeQ.pop();

        // if boundary halfedge, traverse the boundary
        if (isBoundaryHalfedge(pmesh, heh, faceRegionProp)
            && !traversedHalfedges.count(heh))
        {
            std::vector<std::vector<OMPolyMesh::VertexHandle>> foundBoundaries;

            traverseRegionBoundary(pmesh,
                                   heh,
                                   faceRegionProp,
                                   foundBoundaries,
                                   traversedHalfedges);

            for (auto& boundary : foundBoundaries)
            {
                // check if its inner or outer boundary
                if (boundaryArea(boundary, pmesh, glob2plane) > 0)
                {
                    ++outerBoundaryNum;
                    multipolygon.insert(multipolygon.begin(),
                                        std::move(boundary));
                }
                else
                {
                    multipolygon.push_back(std::move(boundary));
                }
            }
        }

        // Check all next halfedges in cw direction
        auto neighboursStart { pmesh.next_halfedge_handle(heh) };
        auto neighbour { neighboursStart };
        do
        {
            if (pmesh.property(faceRegionProp, pmesh.face_handle(neighbour))
                != rIdx)
            {
                break; // exit if not in region
            }
            if (!enqueuedHalfedges.count(neighbour))
            {
                halfedgeQ.push(neighbour);
                enqueuedHalfedges.insert(neighbour);
            }

            neighbour = pmesh.next_halfedge_handle(
                pmesh.opposite_halfedge_handle(neighbour));
        }
        while (neighbour != neighboursStart);
    }

    if (outerBoundaryNum != 1)
    {
        LOGTHROW(err4, std::runtime_error)
            << "Expecting exactly one outer boundary, got: "
            << outerBoundaryNum;
    }
    return multipolygon;
}

/// Create a transformation to 2D coordinates of the face plane
inline math::Matrix4 getGlobal2FacePlane(const OMPolyMesh& pmesh,
                                         const OMPolyMesh::HalfedgeHandle& heh)
{

    auto n { pmesh.calc_face_normal(pmesh.face_handle(heh)) };
    if (!isFinite(n))
    {
        LOGTHROW(err4, std::runtime_error)
            << "Got inf/nan face normal - mesh might contain zero-area faces";
    }
    auto pt { fromOM(pmesh.point(pmesh.from_vertex_handle(heh))) };
    math::Plane3 pl(pt, math::Point3(n[0], n[1], n[2]));
    return math::matrixInvert(math::createPlaneCrs(pl));
}


/// Create multipolygonal faces consisting of connected faces of one region
std::tuple<std::vector<VhMpolyFace>, std::vector<int>>
    getMultiPolygonalFaces(const OMPolyMesh& pmesh,
                           const OMFacePropInt& faceRegionProp)
{
    std::vector<VhMpolyFace> mpolyFaces;
    std::vector<int> mpolyFacesRegions;

    utility::small_set<OMPolyMesh::HalfedgeHandle> traversedHalfedges;
    // check all halfedges
    for (const auto& hehStart : pmesh.halfedges())
    {
        // skip halfedge with no face to the left
        if (!hehStart.is_valid()) { continue; }
        if (hehStart.is_boundary())
        {
            LOG(warn2) << "Mesh is not watertight";
            continue;
        }

        // skip non-boundary
        if (!isBoundaryHalfedge(pmesh, hehStart, faceRegionProp)) { continue; }
        // skip traversed
        if (traversedHalfedges.count(hehStart)) { continue; }

        auto glob2plane { getGlobal2FacePlane(pmesh, hehStart) };

        mpolyFaces.push_back(
            traverseConnectedFacesOfRegion(hehStart,
                                           glob2plane,
                                           pmesh,
                                           faceRegionProp,
                                           traversedHalfedges));
        mpolyFacesRegions.push_back(
            pmesh.property(faceRegionProp, pmesh.face_handle(hehStart)));
    }

    return std::make_tuple(std::move(mpolyFaces), std::move(mpolyFacesRegions));
}

/// Construct multipolygonal mesh based on multipolygonal faces defined by
/// vertex handles in pmesh
MultiPolyMesh<int> createMultiPolyMesh(const OMPolyMesh& pmesh,
                                       const std::vector<VhMpolyFace>& mpFaces,
                                       const std::vector<int>& mpFaceRegions)
{
    MultiPolyMesh<int> omesh;
    std::map<OMPolyMesh::VertexHandle, std::size_t> vertexMap;

    // add all used vertices
    for (const auto& mpFace : mpFaces)
    {
        for (const auto& poly : mpFace)
        {
            for (const auto& v : poly)
            {
                if (vertexMap.count(v)) { continue; }
                auto vIdx { omesh.vertices.size() };
                omesh.vertices.push_back(fromOM(pmesh.point(v)));
                vertexMap[v] = vIdx;
            }
        }
    }

    // add all faces
    omesh.faces.reserve(mpFaces.size());
    for (size_t i = 0; i < mpFaces.size(); i++)
    {
        auto& mpFace { mpFaces[i] };
        geometry::MultiPolyFace face;
        for (const auto& vhPoly : mpFace)
        {
            std::vector<std::size_t> poly;
            for (const auto& v : vhPoly)
            {
                poly.push_back(vertexMap.at(v));
            }
            face.push_back(std::move(poly));
        }
        omesh.faces.push_back(std::move(face));
        omesh.faceLabels.push_back(mpFaceRegions[i]);
    }

    return omesh;
}


} // namespace

MultiPolyMesh<int> polygonizeMesh(const OMPolyMesh& mesh,
                                  const OMFacePropInt& faceRegions)
{
    auto [mpFaces, mpRegions] = getMultiPolygonalFaces(mesh, faceRegions);
    return createMultiPolyMesh(mesh, mpFaces, mpRegions);
}


MultiPolyMesh<int> polygonizeMesh(const Mesh& mesh,
                                  const std::vector<int>& faceRegions)
{
    if (faceRegions.size() != mesh.faces.size())
    {
        LOGTHROW(err4, std::runtime_error)
            << "Face regions contains " << faceRegions.size()
            << " items while the mesh has " << mesh.faces.size() << " faces.";
    }

    // convert to openmesh & add property
    OMPolyMesh pmesh;
    OMFacePropInt faceRegionProp;
    pmesh.add_property(faceRegionProp);

    std::vector<OMPolyMesh::VertexHandle> pVertices(mesh.vertices.size());
    for (size_t i = 0; i < mesh.vertices.size(); i++)
    {
        auto& p { mesh.vertices[i] };
        pVertices[i] = pmesh.add_vertex(OMPolyMesh::Point(p(0), p(1), p(2)));
    }

    for (size_t i = 0; i < mesh.faces.size(); i++)
    {
        auto a { pVertices[mesh.faces[i].a] };
        auto b { pVertices[mesh.faces[i].b] };
        auto c { pVertices[mesh.faces[i].c] };
        auto fh { pmesh.add_face(a, b, c) };
        if (!fh.is_valid())
        {
            LOGTHROW(err4, std::runtime_error) << "Unable to add face";
        }
        pmesh.property(faceRegionProp, fh) = faceRegions[i];
    }

    return polygonizeMesh(pmesh, faceRegionProp);
}

} // namespace geometry

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
 * @file mesh-polygonization.hpp
 * @author Tomas Novak <tomas.novak@melowntech.com>
 *
 * Convert triangular mesh to polygonal mesh.
 */

#ifndef GEOMETRY_MESH_POLYGONIZATION_HPP_INCLUDED_
#define GEOMETRY_MESH_POLYGONIZATION_HPP_INCLUDED_

#include <stddef.h>
#include <vector>

#include "dbglog/dbglog.hpp"
#include "utility/gccversion.hpp"

#include "mesh.hpp"
#include "multipolymesh.hpp"

#ifdef GEOMETRY_HAS_OPENMESH
#    include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#endif


namespace geometry
{

#ifdef GEOMETRY_HAS_OPENMESH

struct OMPolyMeshTraits : public OpenMesh::DefaultTraits
{
    using Point = OpenMesh::Vec3d; // doubles
    FaceAttributes(
        OpenMesh::Attributes::Color
        | OpenMesh::Attributes::Status); // face colors, allow removing faces
    EdgeAttributes(OpenMesh::Attributes::Status); // allow removing edges
    VertexAttributes(OpenMesh::Attributes::Status);
    HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};

/**
 * OpenMesh mesh with polygonal faces
 */
using OMPolyMesh = OpenMesh::PolyMesh_ArrayKernelT<OMPolyMeshTraits>;
using OMFacePropInt = OpenMesh::FPropHandleT<int>;

MultiPolyMesh<int> polygonizeMesh(const OMPolyMesh& mesh,
                                  const OMFacePropInt& faceRegions);

#endif // #ifdef GEOMETRY_HAS_OPENMESH


/**
 * Merge regions of triangles to polygons (with holes).
 *
 * Works on watertight 2-manifolds. Mesh must have correct topology. On errors,
 * check for zero-area faces, non-manifold edges, duplicite vertices, ...
 *
 * NB: OpenMesh is required
 *
 * @param[in] mesh
 * @param[in] faceRegions region index for each face
 * @returns resulting multipolygonal mesh
 */
MultiPolyMesh<int> polygonizeMesh(const Mesh& mesh,
                                  const std::vector<int>& faceRegions)
#ifndef GEOMETRY_HAS_OPENMESH
        UTILITY_FUNCTION_ERROR(
            "Mesh polygonization is available only when compiled with "
            "OpenMesh.")
#endif
    ;

} // namespace geometry


#endif /* GEOMETRY_MESH_POLYGONIZATION_HPP_INCLUDED_ */

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
#ifndef geometry_binmesh_hpp_included_
#define geometry_binmesh_hpp_included_

#include <iosfwd>

#include <boost/filesystem/path.hpp>

#include "./parse-obj.hpp"

namespace geometry {

//! Writes an OBJ mesh in a compact binary format.
void writeBinaryMesh(const boost::filesystem::path &path,
                     const geometry::Obj &mesh);

void writeBinaryMesh(std::ostream &out, const geometry::Obj &mesh);

struct MeshInfo {
    math::Extents3 bbox;
    std::size_t vertexCount;
    std::size_t faceCount;
    std::size_t texCoordCount;
};

//! Loads an OBJ mesh from a compact binary format.
geometry::Obj loadBinaryMesh(const boost::filesystem::path &path
                             , MeshInfo *meshInfo = nullptr);


/** Loads an OBJ mesh from a compact binary format
 * \param path used only for logging
 */
geometry::Obj loadBinaryMesh(std::istream &in
                             , const boost::filesystem::path &path
                             = "unknown"
                             , MeshInfo *meshInfo = nullptr);

} // namespace geometry

#endif // geometry_binmesh_hpp_included_

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

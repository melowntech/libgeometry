/**
 * @file math.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * 3D mesh operations
 */

#ifndef geometry_meshop_hpp_included_
#define geometry_meshop_hpp_included_

#include "utility/gccversion.hpp"

#include "./mesh.hpp"
#include "./parse-obj.hpp"

namespace geometry {

Mesh::pointer simplify(const Mesh &mesh, int faceCount)
#ifndef GEOMETRY_HAS_OPENMESH
    UTILITY_FUNCTION_ERROR("Mesh simplify is available only when compiled with OpenMesh.")
#endif
    ;

Mesh::pointer simplify(const Mesh::pointer &mesh, int faceCount);

/** TODO: remove this once geometry::Obj is no longer used for modeling.
*/
Obj asObj(const Mesh &mesh);

/** TODO: remove this once geometry::Obj is no longer used for modeling.
*/
Obj asObj(const Mesh::pointer &mesh);

void saveAsObj(const Mesh &mesh
               , const boost::filesystem::path &filepath
               , const std::string &mtlName);

void saveAsObj(const Mesh::pointer &mesh
               , const boost::filesystem::path &filepath
               , const std::string &mtlName);

// inline stuff

inline Mesh::pointer simplify(const Mesh::pointer &mesh, int faceCount)
{
    return simplify(*mesh, faceCount);
}

inline Obj asObj(const Mesh::pointer &mesh)
{
    return asObj(*mesh);
}

inline void saveAsObj(const Mesh::pointer &mesh
                      , const boost::filesystem::path &filepath
                      , const std::string &mtlName)
{
    return saveAsObj(*mesh, filepath, mtlName);
}

} // namespace geometry
#endif // geometry_meshop_hpp_included_

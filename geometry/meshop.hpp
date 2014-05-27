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
    UTILITY_FUNCTION_ERROR("Mesh simplification is available only when compiled with OpenMesh.")
#endif
    ;

Mesh::pointer simplify(const Mesh::pointer &mesh, int faceCount);

/** Function that tells how many faces should given cell have.
 */
typedef std::function<std::size_t (const math::Extents2&)> FacesPerCell;

/** Simplify mesh with custom number of faces in cell.
 *
 * \param mesh mesh to simplify
 * \param alignment corner of arbitrary cell in cell grid
 *                  (in mesh-local coordinates)
 * \param cellSize size of cell in cell grid
 * \param facesPerCell assigns face count to cells
 * \return simplified mesh
 */
Mesh::pointer simplifyInGrid(const Mesh &mesh, const math::Point2 &alignment
                             , double cellSize
                             , const FacesPerCell &facesPerCell)
#ifndef GEOMETRY_HAS_OPENMESH
    UTILITY_FUNCTION_ERROR("Mesh simplification is available only when compiled with OpenMesh.")
#endif
    ;

Mesh::pointer simplifyInGrid(const Mesh::pointer &mesh
                             , const math::Point2 &alignment
                             , double cellSize
                             , const FacesPerCell &facesPerCell);

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

void saveAsPly( const Mesh::pointer &mesh
              , const boost::filesystem::path &filepath);

Mesh loadPly(const boost::filesystem::path &filepath);

// inline stuff

inline Mesh::pointer simplify(const Mesh::pointer &mesh, int faceCount)
{
    return simplify(*mesh, faceCount);
}

inline Mesh::pointer simplifyInGrid(const Mesh::pointer &mesh
                                    , const math::Point2 &alignment
                                    , double cellSize
                                    , const FacesPerCell &facesPerCell)
{
    return simplifyInGrid(*mesh, alignment, cellSize, facesPerCell);
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

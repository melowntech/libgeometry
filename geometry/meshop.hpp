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


/** Simplify mesh to certain amount of faces
 * \param mesh mesh to simplify
 * \param faceCount target face count of simplified mesh
 * \return simplified mesh
 */
Mesh::pointer simplify(const Mesh::pointer &mesh, int faceCount);

/** Simplify mesh with maximal geometric error
 * \param mesh mesh to simplify
 * \param maxErr maximal geometric error
 * \return simplified mesh
 */
Mesh::pointer simplifyToError(const Mesh &mesh, double maxErr);

/** Refines mesh. Longest edges are splitted until certain amount of faces is reached
 *
 * \param mesh mesh to refine
 * \param maxFacesCount target faces count of the refined mesh
 * \return refined mesh
 */
Mesh::pointer refine( const Mesh &mesh, uint maxFacesCount);

/** Removes non manifold edges (edges with more than 2 incident faces)
 ** and their incident faces.
 *
 * \param mesh mesh to process
 * \return processed mesh
 */
Mesh::pointer removeNonManifoldEdges( const Mesh& omesh );

/** Clips mesh to the given 3d extents 
 *
 * \param mesh mesh to clip
 * \param extents extents defining where to keep geometry
 * \return processed mesh, texture coordinates are at the moment discarted
 */
Mesh::pointer clip( const Mesh& omesh, const math::Extents3& extents);

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

/** TODO: remove this once geometry::Obj is no longer used for modeling.
*/
Mesh::pointer asMesh(const Obj &obj);



void saveAsObj(const Mesh &mesh
               , const boost::filesystem::path &filepath
               , const std::string &mtlName);

void saveAsObj(const Mesh::pointer &mesh
               , const boost::filesystem::path &filepath
               , const std::string &mtlName);

void saveAsPly( const Mesh::pointer &mesh
              , const boost::filesystem::path &filepath);

void saveAsPly( const Mesh &mesh
               , const boost::filesystem::path &filepath);

Mesh loadPly( const boost::filesystem::path &filepath );

Mesh loadObj( const boost::filesystem::path &filepath );

// inline stuff

inline Mesh::pointer simplify(const Mesh::pointer &mesh, int faceCount)
{
    return simplify(*mesh, faceCount);
}

inline Mesh::pointer simplifyToError(const Mesh::pointer &mesh, double maxErr)
{
    return simplify(*mesh, maxErr);
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

inline void saveAsPly( const Mesh::pointer &mesh
                      , const boost::filesystem::path &filepath){
    return saveAsPly(*mesh, filepath);
}

} // namespace geometry
#endif // geometry_meshop_hpp_included_

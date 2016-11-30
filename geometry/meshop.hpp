/**
 * @file math.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * 3D mesh operations
 */

#ifndef geometry_meshop_hpp_included_
#define geometry_meshop_hpp_included_

#include <boost/optional.hpp>

#include "utility/gccversion.hpp"

#include "./mesh.hpp"
#include "./parse-obj.hpp"

namespace geometry {

enum SimplifyOption
{
    NONE = 0,
    CORNERS =1,
    INNERBORDER = 2,
    OUTERBORDER = 4,
    RMNONMANIFOLDEDGES = 8,
    PREVENTFACEFLIP = 16
};

class SimplifyOptions {
public:
    SimplifyOptions(long flags = 0)
        : flags_(flags), alternativeVertices_() {}

    SimplifyOptions& flags(long flags) { flags_ = flags; return *this; }
    bool check(long flag) const { return flags_ & flag; }

    SimplifyOptions& maxError(const boost::optional<double> &maxError)
    {
        maxError_ = maxError; return *this;
    }
    const boost::optional<double>& maxError() const {
        return maxError_;
    }

    SimplifyOptions& maxEdgeLength(const boost::optional<float> &maxEdgeLength)
    {
        maxEdgeLength_ = maxEdgeLength; return *this;
    }
    const boost::optional<float>& maxEdgeLength() const {
        return maxEdgeLength_;
    }

    SimplifyOptions&
    minAspectRatio(const boost::optional<float> &minAspectRatio)
    {
        minAspectRatio_ = minAspectRatio; return *this;
    }
    const boost::optional<float>& minAspectRatio() const {
        return minAspectRatio_;
    }

    SimplifyOptions&
    alternativeVertices(const math::Points3 *alternativeVertices)
    {
        alternativeVertices_ = alternativeVertices; return *this;
    }
    const math::Points3* alternativeVertices() const {
        return alternativeVertices_;
    }

private:
    long flags_;
    boost::optional<double> maxError_;
    boost::optional<float> maxEdgeLength_;
    boost::optional<float> minAspectRatio_;
    const math::Points3 *alternativeVertices_;
};

Mesh::pointer simplify(const Mesh &mesh, int faceCount
                      , const SimplifyOptions &simplifyOptions
                       =  SimplifyOption::CORNERS
                       | SimplifyOption::RMNONMANIFOLDEDGES)
#ifndef GEOMETRY_HAS_OPENMESH
    UTILITY_FUNCTION_ERROR("Mesh simplification is available only when compiled with OpenMesh.")
#endif
    ;

/** Simplify mesh to certain amount of faces
 * \param mesh mesh to simplify
 * \param faceCount target face count of simplified mesh
 * \return simplified mesh
 */
Mesh::pointer simplify( const Mesh::pointer &mesh, int faceCount
                      , const SimplifyOptions &simplifyOptions
                        = SimplifyOption::CORNERS
                        | SimplifyOption::RMNONMANIFOLDEDGES);

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

/** Removes isolated vertices, e.g vertices incidental with 0 faces 
 * Works with untextured meshes
 *
 * \param mesh mesh to process
 * \return processed mesh
 */
Mesh::pointer removeIsolatedVertices( const Mesh& omesh );

/** Clips mesh to the given 3d extents 
 *
 * \param mesh mesh to clip
 * \param extents extents defining where to keep geometry
 * \return processed mesh, texture coordinates are at the moment discarted
 */
Mesh::pointer clip( const Mesh& omesh, const math::Extents3& extents);

/** Support structure for faces-per-tile calculation.
 */
class FacesPerCell {
public:
    typedef std::function<void (FacesPerCell&)> Functor;

    typedef std::vector<std::size_t> PerCellCount;

    FacesPerCell(const math::Size2_<long> &gridSize, const math::Point2 &origin
                 , const math::Size2f &cellSize)
        : gridSize_(gridSize), origin_(origin), cellSize_(cellSize)
        , goal_(math::area(gridSize))
    {}

    /** Returns cell index.
     */
    std::size_t cellIndex(double x, double y) const;

    std::size_t& operator[](std::size_t index) { return goal_[index]; }
    std::size_t operator[](std::size_t index) const { return goal_[index]; }

    const PerCellCount& get() const { return goal_; }
    std::size_t size() const { return goal_.size(); }

    const math::Size2f& cellSize() const { return cellSize_; }
    math::Extents2 cellExtents(std::size_t index) const;

private:
    const math::Size2_<long> &gridSize_;
    const math::Point2 &origin_;
    const math::Size2f &cellSize_;
    PerCellCount goal_;
};

/** Function that tells how many faces should given cell have.
 */

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
                             , const math::Size2f &cellSize
                             , const FacesPerCell::Functor &facesPerCell
                             , const SimplifyOptions &simplifyOptions
                             =  SimplifyOption::CORNERS
                             | SimplifyOption::RMNONMANIFOLDEDGES)
#ifndef GEOMETRY_HAS_OPENMESH
    UTILITY_FUNCTION_ERROR("Mesh simplification is available only when compiled with OpenMesh.")
#endif
    ;

Mesh::pointer simplifyInGrid(const Mesh::pointer &mesh
                             , const math::Point2 &alignment
                             , const math::Size2f &cellSize
                             , const FacesPerCell::Functor &facesPerCell
                             , const SimplifyOptions &simplifyOptions
                             = SimplifyOption::CORNERS
                             | SimplifyOption::RMNONMANIFOLDEDGES);


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

inline Mesh::pointer simplify(const Mesh::pointer &mesh, int faceCount
                              , const SimplifyOptions &simplifyOptions)
{
    return simplify(*mesh, faceCount, simplifyOptions);
}

inline Mesh::pointer simplifyToError(const Mesh::pointer &mesh, double maxErr)
{
    return simplify(*mesh, maxErr);
}

inline Mesh::pointer simplifyInGrid(const Mesh::pointer &mesh
                                    , const math::Point2 &alignment
                                    , const math::Size2f &cellSize
                                    , const FacesPerCell::Functor &facesPerCell
                                    , const SimplifyOptions &simplifyOptions)
{
    return simplifyInGrid(*mesh, alignment, cellSize, facesPerCell
                          , simplifyOptions);
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

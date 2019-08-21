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
/**
 * @file math.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * 3D mesh operations
 */

#ifndef geometry_meshop_hpp_included_
#define geometry_meshop_hpp_included_

#include <iostream>

#include <boost/optional.hpp>

#include "utility/gccversion.hpp"

#include "mesh.hpp"
#include "parse-obj.hpp"

// Needed to get OpenMesh version
#ifdef GEOMETRY_HAS_OPENMESH
#  include <OpenMesh/Core/System/config.h>
#endif

namespace geometry {

struct GridCell {
    math::Extents2 extent;
    long maxFaceCount;

    GridCell (const math::Extents2 &cellExtent) : extent(cellExtent) {
        maxFaceCount = 0;
    };
};

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

    SimplifyOptions&
    classBasedPlanarisation(const boost::optional<bool> &value)
    {
        classBasedPlanarisation_ = value; return *this;
    }
    const boost::optional<bool>& classBasedPlanarisation() const {
        return classBasedPlanarisation_;
    }

    SimplifyOptions&
    concaveVertexModifier(const boost::optional<float> &value)
    {
        concaveVertexModifier_ = value; return *this;
    }
    const boost::optional<float>& concaveVertexModifier() const {
        return concaveVertexModifier_;
    }

private:
    long flags_;
    boost::optional<double> maxError_;
    boost::optional<float> maxEdgeLength_;
    boost::optional<float> minAspectRatio_;
    const math::Points3 *alternativeVertices_;

    // This is a priority modifier for concave vertexes.
    // Lets denote it K: new_error = old_error * K.
    // For convex and undefined verteces the new_error = old_error.
    // The decimation framework would select verteces
    // with the lowest value of error.
    // Therefore, if K is in (0, 1) then we would prefer
    // these concave verteces (for simplification) over
    // verteces which are not concave.
    // If K in (1, +infinity) then we would have the opposite effect.
    boost::optional<float> concaveVertexModifier_;

    // Use class based planarisation in simplifyToError
    boost::optional<bool> classBasedPlanarisation_;
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
                        | SimplifyOption::RMNONMANIFOLDEDGES)
#ifndef GEOMETRY_HAS_OPENMESH
    UTILITY_FUNCTION_ERROR("Mesh simplification is available only when compiled with OpenMesh.")
#endif
    ;

/** Simplify mesh with maximal geometric error
 * \param mesh mesh to simplify
 * \param maxErr maximal geometric error
 * \return simplified mesh
 */
Mesh::pointer simplifyToError(const Mesh &mesh, double maxErr
                              , const SimplifyOptions &simplifyOptions)
#if !defined(GEOMETRY_HAS_OPENMESH) || (OM_VERSION < 0x60000)
    UTILITY_FUNCTION_ERROR("Error-based mesh simplification is available only when compiled with OpenMesh>=6.")
#endif
    ;

/** Refines mesh. Longest edges are splitted until certain amount of faces is reached
 *
 * \param mesh mesh to refine
 * \param maxFacesCount target faces count of the refined mesh
 * \return refined mesh
 */
Mesh::pointer refine( const Mesh &mesh, unsigned int maxFacesCount);

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
Mesh::pointer clip(const Mesh &omesh, const math::Extents3 &extents);

/** Clips mesh to the given 2d extents 
 *
 * \param mesh mesh to clip
 * \param extents extents defining where to keep geometry
 * \return processed mesh, texture coordinates are at the moment discarted
 */
Mesh::pointer clip(const Mesh &omesh, const math::Extents2 &extents);

/** Appends one mesh to the another. Fixed face indices.
 */
void append(Mesh &mesh, const Mesh &added);

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

void make_gts_class_system_threadsafe()
#ifndef GEOMETRY_HAS_GTS
    UTILITY_FUNCTION_ERROR("Volume-based mesh simplification is available only when compiled with GTS.")
#endif
    ;

Mesh::pointer simplify_gts(const geometry::Mesh &mesh, long edgeCountMax)
#ifndef GEOMETRY_HAS_GTS
    UTILITY_FUNCTION_ERROR("Volume-based mesh simplification is available only when compiled with GTS.")
#endif
    ;

Mesh::pointer simplify_gts_in_grid(const geometry::Mesh &mesh
    , std::vector<std::vector <geometry::GridCell>> &gridCells
    , bool inParallel
    , std::function<math::Point2_<long>(double x, double y)> getGridCell)
#ifndef GEOMETRY_HAS_GTS
    UTILITY_FUNCTION_ERROR("Volume-based mesh simplification is available only when compiled with GTS.")
#endif
    ;

math::Extents2 gridExtents(const math::Extents2 &extents
                           , const math::Point2 &alignment
                           , const math::Size2f &cellSize);
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

/**
 */
struct ObjMaterial {
    std::vector<std::string> libs;
    std::vector<std::string> names;

    ObjMaterial() = default;
    ObjMaterial(const std::string &lib) : libs{lib} {}

    std::string name(std::size_t index) const;
};

void saveAsObj(const Mesh &mesh
               , const boost::filesystem::path &filepath
               , const ObjMaterial &mtl);

void saveAsObj(const Mesh::pointer &mesh
               , const boost::filesystem::path &filepath
               , const ObjMaterial &mtl);

void saveAsObj(const Mesh &mesh, std::ostream &os
               , const ObjMaterial &mtl
               , const boost::filesystem::path &filepath = "UNKNOWN");

void saveAsObj(const Mesh::pointer &mesh, std::ostream &os
               , const ObjMaterial &mtl
               , const boost::filesystem::path &filepath = "UNKNOWN");

void saveAsPly( const Mesh::pointer &mesh
              , const boost::filesystem::path &filepath);

void saveAsPly( const Mesh &mesh
               , const boost::filesystem::path &filepath);

Mesh loadPly( const boost::filesystem::path &filepath );

Mesh loadObj(const boost::filesystem::path &filepath
             , ObjMaterial *mtl = nullptr);

// parses PLY from a file, throws on any error
void loadPly(ObjParserBase &parser, const boost::filesystem::path &path);

struct MeshInfo {
    std::size_t vertexCount;
    std::size_t faceCount;

    MeshInfo(std::size_t vertexCount, std::size_t faceCount)
        : vertexCount(vertexCount), faceCount(faceCount)
    {}
};

MeshInfo measurePly(const boost::filesystem::path &path);

/** Splits mesh by image identifier. Duplicates vertices shared by different
 *  images.
 */
Mesh::list splitById(const Mesh &mesh);

// inline stuff

inline Mesh::pointer simplify(const Mesh::pointer &mesh, int faceCount
                              , const SimplifyOptions &simplifyOptions)
{
    return simplify(*mesh, faceCount, simplifyOptions);
}

inline Mesh::pointer simplifyToError(const Mesh::pointer &mesh, double maxErr, const SimplifyOptions &options)
{
    return simplifyToError(*mesh, maxErr, options);
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
                      , const ObjMaterial &mtl)
{
    return saveAsObj(*mesh, filepath, mtl);
}

inline void saveAsObj(const Mesh::pointer &mesh, std::ostream &os
                      , const ObjMaterial &mtl
                      , const boost::filesystem::path &filepath)
{
    return saveAsObj(*mesh, os, mtl, filepath);
}

inline void saveAsPly( const Mesh::pointer &mesh
                      , const boost::filesystem::path &filepath)
{
    return saveAsPly(*mesh, filepath);
}

} // namespace geometry
#endif // geometry_meshop_hpp_included_

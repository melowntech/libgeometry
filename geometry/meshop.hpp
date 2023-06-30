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
#include <boost/container/small_vector.hpp>

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
    long long maxFaceCount;

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
    SimplifyOptions(long long flags = 0)
        : flags_(flags), alternativeVertices_() {}

    SimplifyOptions& flags(long long flags) { flags_ = flags; return *this; }
    bool check(long long flag) const { return flags_ & flag; }

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

    SimplifyOptions& vertexWeights(const boost::optional<std::vector<float>>& weights) {
        weights_ = weights;
        return *this;
    }
    const boost::optional<std::vector<float>>& vertexWeights() const {
        return weights_;
    }

private:
    long long flags_;
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

    // Multiplier of quadric cost for each vertex of the mesh
    boost::optional<std::vector<float>> weights_;
};


template <typename T>
struct EdgeKeyTmpl
{
    T v1_, v2_; // vertex indices

    EdgeKeyTmpl(T v1, T v2)
    {
        v1_ = std::min(v1, v2);
        v2_ = std::max(v1, v2);
    }

    bool operator<(const EdgeKeyTmpl& other) const
    {
        return (v1_ == other.v1_) ? (v2_ < other.v2_) : (v1_ < other.v1_);
    }
};


using EdgeKey = EdgeKeyTmpl<Face::index_type>;
using NonManifoldEdge = boost::container::small_vector<Face::index_type, 2>;
using EdgeMap = std::map<EdgeKey, NonManifoldEdge>;


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

/**
 * @brief Simplifies the mesh in place to certain amount of faces
 * @param mesh mesh to simplify
 * @param faceCount target face count
 * @param simplifyOptions additional parameters of the simplification
 */
void simplifyInPlace(Mesh &mesh, int faceCount
                    , const SimplifyOptions &simplifyOptions
                     =  SimplifyOption::CORNERS
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


EdgeMap getNonManifoldEdgeMap(const Mesh& mesh);


/** Removes non manifold edges (edges with more than 2 incident faces)
 ** and their incident faces.
 *
 * \param mesh mesh to process
 * \return processed mesh
 */
Mesh::pointer removeNonManifoldEdges(Mesh omesh);

/**
 * Face incidency - for each face lists all neighbors, supports non-manifolds
 */
using FaceFaceTable = std::vector<boost::container::small_vector<std::size_t, 3>>;

/** Computes face incidency map, supports non manifold edges (face may have >3
 * neighbors).
 *
 * NB: when two faces are incident over multiple edges, they are repeated
 * respective amount of times
 *
 * @param[in] mesh mesh to process
 * @returns face incidency table
 */
FaceFaceTable getFaceFaceTableNonManifold(const Mesh& mesh);

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
Mesh clip(const Mesh &omesh, const math::Extents3 &extents);

/** Clips mesh to the given 2d extents
 *
 * \param mesh mesh to clip
 * \param extents extents defining where to keep geometry
 * \return processed mesh, texture coordinates are at the moment discarted
 */
Mesh clip(const Mesh &omesh, const math::Extents2 &extents);

/** Clips mesh to the given 2d contour
 *
 * \param mesh mesh to clip
 * \param clipRegion clipping region
 * \return processed mesh, including texture coordinates
 */
Mesh clip(const Mesh &omesh, const math::MultiPolygon &clipRegion);

/** Appends one mesh to the another. Fixed face indices.
 */
void append(Mesh &mesh, const Mesh &added);

/** Support structure for faces-per-tile calculation.
 */
class FacesPerCell {
public:
    typedef std::function<void (FacesPerCell&)> Functor;

    typedef std::vector<std::size_t> PerCellCount;

    FacesPerCell(const math::Size2ll &gridSize
                 , const math::Point2 &origin
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
    const math::Size2ll &gridSize_;
    const math::Point2 &origin_;
    const math::Size2f &cellSize_;
    PerCellCount goal_;
};

Mesh::pointer simplify_gts(const geometry::Mesh &mesh, long long edgeCountMax
                           , double costMax)
#ifndef GEOMETRY_HAS_GTS
    UTILITY_FUNCTION_ERROR("Volume-based mesh simplification is available only when compiled with GTS.")
#endif
    ;

Mesh::pointer simplify_gts_in_grid(const geometry::Mesh &mesh
    , std::vector<std::vector <geometry::GridCell>> &gridCells
    , bool inParallel
    , std::function<math::Point2ll(double x, double y)> getGridCell)
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
    ObjMaterial(const char *lib) : libs{lib} {}

    std::string name(std::size_t index) const;
};

/** Callback for setting up stream before saving mesh as OBJ.
 *  Returns if stream precision has been set.
 *  Stream precision is configured when false is returned.
 */
struct ObjStreamSetup {
    virtual bool vertex(std::ostream&) const { return false; };
    virtual bool tx(std::ostream&) const { return false; };

    virtual ~ObjStreamSetup() = default;
};

void saveAsObj(const Mesh &mesh
               , const boost::filesystem::path &filepath
               , const ObjMaterial &mtl
               , const ObjStreamSetup &streamSetup = ObjStreamSetup());

void saveAsObj(const Mesh &mesh
               , const boost::filesystem::path &filepath
               , const ObjMaterial &mtl
               , const std::function<bool(std::ostream&)> &streamSetup);


void saveAsObj(const Mesh::pointer &mesh
               , const boost::filesystem::path &filepath
               , const ObjMaterial &mtl);

void saveAsObj(const Mesh &mesh, std::ostream &os
               , const ObjMaterial &mtl
               , const boost::filesystem::path &filepath = "UNKNOWN"
               , bool setFormat = true);

void saveAsObj(const Mesh::pointer &mesh, std::ostream &os
               , const ObjMaterial &mtl
               , const boost::filesystem::path &filepath = "UNKNOWN"
               , bool setFormat = true);

void saveAsObj(const Mesh &mesh, std::ostream &out
               , const ObjMaterial &mtl
               , const boost::filesystem::path &filepath
               , const ObjStreamSetup &streamSetup);

void saveAsGzippedObj(const Mesh &mesh
                      , const boost::filesystem::path &filepath
                      , const ObjMaterial &mtl
                      , const ObjStreamSetup &streamSetup = ObjStreamSetup())
#ifndef GEOMETRY_HAS_BIO
    UTILITY_FUNCTION_ERROR("Gzip support is available only when compiled with Boost.IOStreams.")
#endif
    ;

void saveAsGzippedObj(const Mesh &mesh
                      , const boost::filesystem::path &filepath
                      , const ObjMaterial &mtl
                      , const std::function<bool(std::ostream&)>
                      &streamSetup = {})
#ifndef GEOMETRY_HAS_BIO
    UTILITY_FUNCTION_ERROR("Gzip support is available only when compiled with Boost.IOStreams.")
#endif
    ;

void saveAsPly( const Mesh::pointer &mesh
              , const boost::filesystem::path &filepath);

void saveAsPly( const Mesh &mesh
               , const boost::filesystem::path &filepath);

Mesh loadPly( const boost::filesystem::path &filepath );

// parses PLY from a file, throws on any error
void loadPly(ObjParserBase &parser, const boost::filesystem::path &path);

Mesh loadObj(const boost::filesystem::path &filepath
             , ObjMaterial *mtl = nullptr);

void loadObj(ObjParserBase &parser, const boost::filesystem::path &filepath);

struct MeshInfo {
    std::size_t vertexCount;
    std::size_t faceCount;

    MeshInfo(std::size_t vertexCount, std::size_t faceCount)
        : vertexCount(vertexCount), faceCount(faceCount)
    {}
};

MeshInfo measurePly(const boost::filesystem::path &path);

bool objHasVertexOrFace(const boost::filesystem::path &path);

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
                      , const boost::filesystem::path &filepath
                      , bool setFormat)
{
    return saveAsObj(*mesh, os, mtl, filepath, setFormat);
}

inline void saveAsPly( const Mesh::pointer &mesh
                      , const boost::filesystem::path &filepath)
{
    return saveAsPly(*mesh, filepath);
}

namespace detail {

class CallbackStreamSetup : public ObjStreamSetup {
public:
    CallbackStreamSetup(const std::function<bool(std::ostream&)> &func)
        : func_(func)
    {}

    virtual bool vertex(std::ostream &os) const {
        return func_(os);
    }

    virtual bool tx(std::ostream &os) const {
        return func_(os);
    }

private:
    std::function<bool(std::ostream&)> func_;
};

} // namespace detail

inline void saveAsObj(const Mesh &mesh
                      , const boost::filesystem::path &filepath
                      , const ObjMaterial &mtl
                      , const std::function<bool(std::ostream&)> &streamSetup)
{
    saveAsObj(mesh, filepath, mtl
              , detail::CallbackStreamSetup(streamSetup));
}

inline void saveAsGzippedObj(const Mesh &mesh
                             , const boost::filesystem::path &filepath
                             , const ObjMaterial &mtl
                             , const std::function<bool(std::ostream&)>
                             &streamSetup)
{
    saveAsGzippedObj(mesh, filepath, mtl
                     , detail::CallbackStreamSetup(streamSetup));
}

} // namespace geometry
#endif // geometry_meshop_hpp_included_

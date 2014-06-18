/**
 * @file mesh-voxelizer.hpp
 * @author Tomas Drinovsky <ondrej.prochazka@citationtech.net>
 *
 * Mesh voxelization class.
 *
 * This module provides a class voxelization of the generic meshes.
 * The implementation is based on the paper "Simplification and Repair of
 * Polygonal Models Using Volumetric Techniques" by F.S. Nooruddin and Greg Turk
 */



#ifndef GEOMETRY_MESH_VOXELIZER_HPP
#define GEOMETRY_MESH_VOXELIZER_HPP

#include "geometry/mesh.hpp"
#include "geometry/volume.hpp"

#include <memory>

#include <boost/filesystem.hpp>

namespace geometry{

namespace fs = boost::filesystem;

/**
 * @brief MeshVoxelizer
 * @details Class Mesh Voxelizer is able to voxelize arbitrary mesh
 * (even non-water-tight), refilter the voxel grid to different resolution
 * and extract iso surface of this voxel grid.
 * Right now the voxel sizes must be all same - voxel must be a cube.
 */

class MeshVoxelizer{

public:
    enum class Method{ PARITY_COUNT, RAY_STABING };

    struct Parameters {
        float voxelSize;
        bool addFloor;
        float isoThreshold;
        float sealFactor;
        Method method;

        Parameters():
            voxelSize(0.25)
          , addFloor(true)
          , isoThreshold(0.5)
          , sealFactor(3)
          , method(Method::PARITY_COUNT)
         {};
    };

    /**
     * @brief Constructs MeshVoxelizer class
     * @param mesh The geometry to voxelize
     * @param params parameters of MeshVoxelizer
     */

    MeshVoxelizer(const Parameters & params);

    void add( geometry::Mesh & mesh );
    void voxelize();
    std::shared_ptr<ScalarField_t<float>> volume();
    void reset();


private:
    struct LayeredZBuffer{
        typedef std::vector<std::vector<std::list<float>>> LZBData;
        math::Size2i size;
        LZBData data;

        LayeredZBuffer(math::Size2i & size){
            this->size = size;
            data = LZBData(size.width);
            for(int i=0;i<size.width; ++i){
                data[i]= std::vector<std::list<float>>(size.height);
            };
        };
    };

    struct ProjectionResult{
        math::Matrix4 transformation;
        LayeredZBuffer buffer;

        ProjectionResult(math::Matrix4 transformation, LayeredZBuffer buffer):
            transformation(transformation), buffer(buffer){};
    };


    typedef std::vector<ProjectionResult> ProjectionResults;

    Parameters params_;
    std::shared_ptr<ScalarField_t<float>> volume_;
    std::vector<geometry::Mesh> meshes;
    ProjectionResults projections;

    math::Matrix4 orthoProjMat( const math::Point3 &direction
                               , const math::Extents3 &extents
                               , const float &voxelSize
                               , math::Size2 &viewport);

    void addSealToMesh(geometry::Mesh & mesh);
    void fillVolumeFromSeal();

    /**
     * @brief rasterizeMesh renders the mesh into layered z-buffer using
     * given projection matrix.
     * @param mesh
     * @param projection
     */
    void rasterizeMesh( const geometry::Mesh & mesh
                       , const math::Matrix4 & projMat
                       , LayeredZBuffer & lZBuffer);

    /**
     * @brief Determines if point on given position is inside.
     * From each projection it determines if the voxel is inside or outside
     * by parity-counting or by ray-stabing method.
     * @return True if the voxel is inside
     */
    bool isInside(  const math::Point3 & position
                  , const ProjectionResults & projectionResults);

    void visualizeDepthMap(const ProjectionResult &proj
                           , const math::Extents3 &extents, const fs::path &path);

};

} //namespace geometry

#endif// GEOMETRY_MESH_VOXELIZER_HPP

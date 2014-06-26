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
    typedef ScalarField_t<unsigned short, VolumeArray<unsigned short>> Volume;


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
    std::shared_ptr<Volume> volume();
    void reset();
private:
    struct LayeredZBuffer{
        typedef std::vector<std::vector<std::vector<float>>> LZBData;
        typedef std::vector<float>::iterator Iterator;
        math::Size2i size;
        LZBData data;

        LayeredZBuffer(math::Size2i & size){
            this->size = size;
            data = LZBData(size.width);
            for(int i=0;i<size.width; ++i){
                data[i]= std::vector<std::vector<float>>(size.height);
                for(int j=0;j<size.height; ++j){
                    data[i][j].reserve(0);
                };
            };
        }

        void sortCells(){
            for( auto &col : data){
                for( auto &cell : col){
                    std::sort(cell.begin(),cell.end());
                }
            }
        }

        void checkCells(){
            for( uint col=0; col < data.size(); col++){
                for( uint cell=0; cell < data[col].size(); cell++){
                    if(data[col][cell].size()>100){
                        std::cout << col<<" : "<<cell
                                  <<" - " << data[col][cell].size()<<std::endl;
                    }
                }
            }
        }

        long avgCellSize(){
            long totalSize=0;
            for( uint col=0; col < data.size(); col++){
                for( uint cell=0; cell < data[col].size(); cell++){
                    totalSize += data[col][cell].size();
                }
            }
            return totalSize/(size.width*size.height);
        }

        long mem(){
            long dataMem=0;
            for( auto &col : data){
                for( auto &cell : col){
                    dataMem+= cell.size()*sizeof(float);
                }
            }

            return sizeof(std::vector<float>)*(size.width*size.height+size.width)
                     + dataMem;
        }


        Iterator begin(uint x, uint y){
            return data[x][y].begin();
        }

        Iterator end(uint x, uint y){
            return data[x][y].end();
        }
    };

    struct CompressedLayeredZBuffer{
        typedef std::vector<float>::iterator Iterator;
        math::Size2i size;
        std::vector<float> data;
        std::vector<size_t> pos;
        std::vector<unsigned short> count;

        CompressedLayeredZBuffer():size(math::Size2(0,0)){}

        CompressedLayeredZBuffer(LayeredZBuffer &lzBuffer){
            this->size = lzBuffer.size;
            pos.reserve(size.width*size.height);
            count.reserve(size.width*size.height);

            std::size_t size = 0;
            for( auto &col : lzBuffer.data){
                for( auto &cell : col){
                    count.push_back(cell.size());
                    pos.push_back(size);
                    size+=cell.size();
                }
            }
            data.reserve(size);
            for( auto &col : lzBuffer.data){
                for( auto &cell : col){
                    for( auto &zval : cell){
                        data.push_back(zval);
                    }
                }
            }
        }

        Iterator begin(uint x, uint y){
            return data.begin()+pos[(std::size_t)x*size.height+y];
        }

        Iterator end(uint x, uint y){
            return data.begin()+pos[(std::size_t)x*size.height+y]
                        +count[(std::size_t)x*size.height+y];
        }

    };

    struct Projection{
        math::Matrix4 transformation;
        math::Size2 viewportSize;

        Projection( math::Matrix4 transformation, math::Size2 viewportSize):
                  transformation(transformation),viewportSize(viewportSize){}
    };

    struct ProjectionResult{
        math::Matrix4 transformation;
        CompressedLayeredZBuffer buffer;

        ProjectionResult( math::Matrix4 transformation
                        , CompressedLayeredZBuffer buffer):
            transformation(transformation), buffer(buffer){};
    };

    typedef std::vector<ProjectionResult> ProjectionResults;
    typedef std::vector<Projection> Projections;
    typedef std::vector<CompressedLayeredZBuffer> CompLZBuffers;

    Parameters params_;
    std::shared_ptr<Volume> volume_;
    std::vector<geometry::Mesh*> meshes;

    Projection orthoProj( const math::Point3 &direction
                               , const math::Extents3 &extents
                               , const float &voxelSize);

    geometry::Mesh sealOfMesh(geometry::Mesh & mesh);
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
                  , ProjectionResults & projectionResults);
/*
    void visualizeDepthMap(const ProjectionResult &proj
                           , const math::Extents3 &extents, const fs::path &path);
*/
};

} //namespace geometry

#endif// GEOMETRY_MESH_VOXELIZER_HPP

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

#include "utility/openmp.hpp"

#include "geometry/meshop.hpp"
#include "geometry/mesh-voxelizer.hpp"
#include "math/geometry_core.hpp"
#include "imgproc/scanconversion.hpp"

//#include <opencv2/opencv.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <algorithm>
//#define RASTERIZE_MESH_DEBUG



namespace geometry{

MeshVoxelizer::MeshVoxelizer(const Parameters &params){
    params_ = params;
}

void MeshVoxelizer::add( geometry::Mesh & mesh ){
    meshes.push_back(&mesh);
}

void MeshVoxelizer::voxelize(){
    if(meshes.size()<=0){
        LOG( warn2 )<<"Zero meshes to voxelize. Skipping voxelization.";
        return;
    }

    //add floor to mimic closed mesh
    std::vector<geometry::Mesh> seals;
    if(params_.addSeal){
        for(auto & mesh: meshes){
            seals.push_back(sealOfMesh(*mesh));
        }
    }

    //compute united extents
    math::Extents3 extents(math::InvalidExtents{});
    for(auto & mesh: meshes){
        math::Extents3 meshExtents = math::computeExtents(mesh->vertices);
        extents = math::unite(extents, meshExtents);
    }
    for(auto & mesh: seals){
        math::Extents3 meshExtents = math::computeExtents(mesh.vertices);
        extents = math::unite(extents, meshExtents);
    }

    LOG( info2 )<<"Meshes extents: "<<extents;

    // volume size

    {
        math::Point3 extentsSize = extents.ur-extents.ll;
        math::Point3i volumeGridRes(
                    std::ceil(extentsSize(0)/params_.voxelSize)
                    ,std::ceil(extentsSize(1)/params_.voxelSize)
                    ,std::ceil(extentsSize(2)/params_.voxelSize));
        math::Point3 newExtentsSize = volumeGridRes*params_.voxelSize;
        math::Point3 extentsCenter = math::center(extents);

        extents = math::Extents3( extentsCenter-newExtentsSize/2
                                 , extentsCenter+newExtentsSize/2);

        if ( valid(params_.overrideExtents) ) {
            LOG(info2) << "Using supplied extents for voxelization.";

            extents = params_.overrideExtents;

            for (uint i(0); i < 3; ++i) {
                volumeGridRes(i)
                    = std::ceil(size(extents)(i)/params_.voxelSize);
            }
        }

        LOG(info2) << "Volumetric grid resolution: " << volumeGridRes;
        LOG(info2) << "Volume extents: " << extents;
    }


    //generate directions for 3 main axes + axes of icosahedron faces
    std::vector<math::Point3> directions;
    directions.push_back(-math::Point3(1,0,0));
    directions.push_back(-math::Point3(0,1,0));
    directions.push_back(-math::Point3(0,0,1));

    directions.push_back(-math::Point3(1,1,1));
    directions.push_back(-math::Point3(1,-1, 1));
    directions.push_back(-math::Point3(-1,1,1));
    directions.push_back(-math::Point3(-1,-1,1));
    directions.push_back(-math::Point3(1,0,1));
    directions.push_back(-math::Point3(-1,0,1));
    directions.push_back(-math::Point3(0,1,1));
    directions.push_back(-math::Point3(0,-1,1));
    directions.push_back(-math::Point3(1,1,0));
    directions.push_back(-math::Point3(1,-1,0));

    //generate projection matrices and layered zbuffers

    Projections projections;
    ProjectionResults results;
    for(const auto & dir: directions){
        Projection projection = orthoProj( dir, extents, params_.voxelSize);
        projections.push_back(projection);
    }

    LOG( info2 )<<"Rasterizing meshes";
    //rasterize mesh using generated projection matrices
    for(uint p=0; p<projections.size(); ++p){
        LayeredZBuffer buffer(projections[p].viewportSize);
        for(auto & mesh: meshes){
            rasterizeMesh(*mesh, projections[p].transformation, buffer);
        }
        for(auto & mesh: seals){
            rasterizeMesh(mesh, projections[p].transformation, buffer);
        }
        buffer.sortCells();
        results.push_back(ProjectionResult( projections[p].transformation
                                      , CompressedLayeredZBuffer(buffer)));
    }

#ifdef RASTERIZE_MESH_DEBUG
    uint i=0;
    for( auto &proj : results){
        fs::path path = fs::path(std::string("depth")
                                 +boost::lexical_cast<std::string>(i++)
                                 +std::string(".png"));
        visualizeDepthMap(proj, extents ,path);
    }
#endif



    //create new volume
    /**************
     * Ugly hack (shaveVolume)
     *
     * Voxelization has problems with edges - create volume 2*cells
     * smaller from both x and y directions
     *
     * Voxelization itself should be fixed instead
     **************************/
    {
        math::Size3i volSize( std::ceil( size(extents)(0)/params_.voxelSize )
                            , std::ceil( size(extents)(1)/params_.voxelSize )
                            , std::ceil( size(extents)(2)/params_.voxelSize ));

        math::Point3 ll( extents.ll(0)
                       , extents.ll(1)
                       , extents.ll(2));

        if (params_.shaveVolume) {
            for (uint i(0); i < 2; ++i) {
                // strip 2 from each side
                volSize(i) -= 4;
                // shift extents's ll by 2 voxel size because that is what we stripped
                ll(i) += 2*params_.voxelSize;
            }
        }

        // prepare volume for filtering and downsampling: reserve space so all
        // dimensions can be inflated to closest larger odd value
        // add few layers of cells on the ceiling of the volume as a margin
        // for subsampling.
        math::Size3i capacity( volSize(0) + !(volSize(0) % 2)
                             , volSize(1) + !(volSize(1) % 2)
                             , volSize(2) + !(volSize(2) % 2) + 4);

        LOG(info2) << "Creating volume of size " << volSize
                   << " and capacity " << capacity;


        volume_ = std::unique_ptr<Volume>( new Volume( ll
                                                     , params_.voxelSize
                                                     , volSize
                                                     , VoxelizerUnit::empty()
                                                     , capacity));
    }

    math::Size3i vSize = volume_->cSize();
    long long volMem = (long long)vSize.width * vSize.height
                    * vSize.depth * sizeof(unsigned short);
    LOG( info2 )<<"Memory consumption of volume: "
        << volMem/1024.0/1024.0/1024.0 << " GB.";

    long long mem = 0;
    for(auto res : results){
        mem += res.buffer.mem();
    }

    LOG( info2 )<<"Memory consumption of zbuffers: "
        << mem/1024.0/1024.0/1024.0 << " GB.";


    uint progress = 0;
    LOG( info2 )<<"Voxelization progress: "<<progress;
#ifdef _OPENMP
    #pragma omp parallel for schedule( dynamic, 10 )
#endif
    for(int x=0; x< vSize.width; ++x){
        for(int y=0; y< vSize.height; ++y){
            for(int z=0; z< vSize.depth; ++z){
                //voxel center position
                auto voxCenter
                        = volume_->grid2geo(math::Point3(x,y,z));
                if(isInside({voxCenter.x, voxCenter.y, voxCenter.z}, results)){
                        volume_->set( x, y, z,VoxelizerUnit::full());
                }

            }
        }
        //voxelization progress
        uint newProgress = std::ceil((float)x/vSize.width * 100);
        if(newProgress>progress){
            progress=newProgress;
            LOG( info2 )<<"Voxelization progress: "<<progress;
        }
    }

    if(params_.addSeal){
        LOG( info2 )<<"Filling volume from seal";
        fillVolumeFromSeal();
    }
}

std::shared_ptr<MeshVoxelizer::Volume> MeshVoxelizer::volume(){
    return volume_;
}

void MeshVoxelizer::reset(){
    std::vector<geometry::Mesh*>().swap(meshes);
}


MeshVoxelizer::Projection MeshVoxelizer::orthoProj(const math::Point3 &direction
                            , const math::Extents3 &extents
                            , const float &voxelSize
                            ){
    LOG( debug )<<"Projection direction: "<<direction;
    LOG( debug )<<"Mesh extents: "<<extents;

    math::Point3 normalizedDir = math::normalize(direction);

    math::Matrix4 projMat = ublas::identity_matrix<double>(4);
    //translate so the center of the extents is in the middle
    math::Point3 extentsCenter = math::center(extents);
    projMat(0,3) = -extentsCenter[0];
    projMat(1,3) = -extentsCenter[1];
    projMat(2,3) = -extentsCenter[2];

    //LOG( debug )<<"Projection matrix - after centering: "<<projMat;

    //rotate scene in specified direction
    math::Matrix4 rotation = ublas::identity_matrix<double>(4);

    math::Point3 up(0,0,-1);
    if(std::abs(ublas::inner_prod(normalizedDir,math::Point3(0,0,1)))==1){
        up = math::Point3(0,1,0);
    }
    math::Point3 right = math::normalize(math::crossProduct(up,direction));
    up = math::normalize(math::crossProduct(normalizedDir,right));

    auto c1 = ublas::row(rotation,0);
    auto c2 = ublas::row(rotation,1);
    auto c3 = ublas::row(rotation,2);

    ublas::subrange(c1,0,3) = right;
    ublas::subrange(c2,0,3) = up;
    ublas::subrange(c3,0,3) = normalizedDir;

    projMat = prod(rotation,projMat);

    //find out extents images
    math::Extents2 viewExtents(math::InvalidExtents{});
    math::Points3 extPoints;
    extPoints.push_back(math::Point3(extents.ll[0],extents.ll[1],extents.ll[2]));
    extPoints.push_back(math::Point3(extents.ll[0],extents.ll[1],extents.ur[2]));
    extPoints.push_back(math::Point3(extents.ll[0],extents.ur[1],extents.ll[2]));
    extPoints.push_back(math::Point3(extents.ll[0],extents.ur[1],extents.ur[2]));
    extPoints.push_back(math::Point3(extents.ur[0],extents.ll[1],extents.ll[2]));
    extPoints.push_back(math::Point3(extents.ur[0],extents.ll[1],extents.ur[2]));
    extPoints.push_back(math::Point3(extents.ur[0],extents.ur[1],extents.ll[2]));
    extPoints.push_back(math::Point3(extents.ur[0],extents.ur[1],extents.ur[2]));

    for(const auto &point : extPoints){
        math::Point3 pointImg = transform(projMat,point);
        viewExtents.ll[0] = std::min(pointImg[0],viewExtents.ll[0]);
        viewExtents.ll[1] = std::min(pointImg[1],viewExtents.ll[1]);
        viewExtents.ur[0] = std::max(pointImg[0],viewExtents.ur[0]);
        viewExtents.ur[1] = std::max(pointImg[1],viewExtents.ur[1]);
    }

    //scale scene so the extents fit to space < -1,1 > CC/NDC
    math::Matrix4 scaleMat = ublas::identity_matrix<double>(4);
    scaleMat(0,0) = 2./size(viewExtents).width;
    scaleMat(1,1) = 2./size(viewExtents).height;

    projMat = prod(scaleMat, projMat);

    //scale and transform scene to VC
    math::Matrix4 transformMat = ublas::identity_matrix<double>(4);
    transformMat(0,3)=1;
    transformMat(1,3)=1;

    //calculate viewport size
    LOG( debug ) << "View extents: "<< viewExtents;

    math::Size2 viewport( std::ceil(size(viewExtents).width/voxelSize)
                        , std::ceil(size(viewExtents).height/voxelSize));

    LOG( debug ) << "Viewport:" << viewport;

    scaleMat = ublas::identity_matrix<double>(4);
    scaleMat(0,0) = viewport.width/2;
    scaleMat(1,1) = viewport.height/2;

    projMat = prod(transformMat, projMat);
    projMat = prod(scaleMat, projMat);

    LOG( debug )<<"Projection matrix - final from WC to CC: "<<projMat;

    for(const auto &point : extPoints){
        math::Point3 pointImg = transform(projMat,point);
        LOG( debug ) << point <<" - "<< pointImg;
    }

    return Projection(projMat, viewport);
}

geometry::Mesh MeshVoxelizer::sealOfMesh(geometry::Mesh & mesh){

    math::Extents3 extents = math::computeExtents(mesh.vertices);

    float offset = params_.voxelSize * params_.sealFactor;

    geometry::Mesh seal;
    math::Point3 extSize = extents.ur-extents.ll;

    double cellWidth = offset;

    int cols=std::ceil(extSize(0)/cellWidth);
    int rows=std::ceil(extSize(1)/cellWidth);

    std::vector<std::vector<double>> minMap
        = std::vector<std::vector<double>>(cols);
    for(int x=0; x<cols; ++x){
        minMap[x] = std::vector<double>(rows,INFINITY);
    }

    for(const auto& vertex: mesh.vertices){
        int x = std::min((int)std::floor((vertex(0)-extents.ll(0))/cellWidth),cols-1);
        int y = std::min((int)std::floor((vertex(1)-extents.ll(1))/cellWidth),rows-1);
        minMap[x][y]=std::min(minMap[x][y],vertex(2));
    }

    for(int x=0; x< cols; ++x){
        for(int y=0; y< rows; ++y){
            if(x==0 || y==0 || x==cols-1 || y==rows-1){
                //fill unsetted values
                int searchOffset=1;
                while(!std::isfinite(minMap[x][y])){
                    for(int xoff=-searchOffset; xoff<searchOffset; ++xoff ){
                        for(int yoff=-searchOffset; yoff<searchOffset; ++yoff ){
                            if( x+xoff>0 && x+xoff<(int)cols
                               && y+yoff>0 && y+yoff<(int)rows){
                                minMap[x][y]=std::min( minMap[x][y]
                                                     , minMap[x+xoff][y+yoff]);
                            }
                        }
                    }
                    searchOffset++;
                    if(searchOffset>params_.sealFactor*2){
                        minMap[x][y] = extents.ll(2);
                    }
                }
            }
        }
    }

    for(int x=1; x< cols-1; ++x){
        for(int y=1; y< rows-1; ++y){
            minMap[x][y] = extents.ll(2);
        }
    }

    for(int y=0; y< rows; ++y){
        for(int x=0; x< cols; ++x){
            seal.vertices.push_back(
                math::Point3( x*cellWidth+extents.ll(0)
                            , y*cellWidth+extents.ll(1)
                            , minMap[x][y]-offset ));
        }
    }

    for(int x=0; x< cols-1; ++x){
        for(int y=0; y< rows-1; ++y){
            seal.addFace(x+(y*cols), x+1+((y+1)*cols), x+1+(y*cols));
            seal.addFace(x+(y*cols), x+((y+1)*cols), x+1+((y+1)*cols));
        }
    }

    return seal;
}

void MeshVoxelizer::fillVolumeFromSeal(){
    math::Size3i vSize = volume_->cSize();
    for(int x = 0;x<vSize.width; ++x){
        for(int y = 0;y<vSize.height; ++y){
            //find first full voxel from the direction of seal
            int full = -1;
            for(int z = 0;z<vSize.depth; ++z){
                if(volume_->get(x,y,z)!=VoxelizerUnit::empty()){
                        full = z;
                        break;
                }
            }
            //fill from bottom up
            for(int z = 0;z<full; ++z){
                volume_->set(x,y,z,VoxelizerUnit::full());
            }
        }
    }
}

void MeshVoxelizer::rasterizeMesh( const Mesh &mesh
                                 , const math::Matrix4 &projMat
                                 , LayeredZBuffer & lZBuffer){
    std::vector<imgproc::Scanline> scanlines;

    LOG(info2) << "rasterizing " << mesh.faces.size() << " triangles";

    // draw all faces into the zBuffer
    for (const auto &face : mesh.faces)
    {
        cv::Point3f tri[3];

        int i(0);
        for (auto it(mesh.begin(face)); it != mesh.end(face); ++it, ++i) {
            math::Point3d pt(transform(projMat, *it));
            tri[i] = {float(pt(0)), float(pt(1)), float(pt(2))};
        }

        scanlines.clear();
        imgproc::scanConvertTriangle(tri, 0, lZBuffer.size.height, scanlines);

        for (const auto& sl : scanlines)
        {
            imgproc::processScanline(sl, 0, lZBuffer.size.width,
                [&](int x, int y, float z) {
                    lZBuffer.data[x][y].push_back(z);
                } );
        }
    }


}

bool MeshVoxelizer::isInside( const math::Point3 & position
                      ,  ProjectionResults & projectionResults){
    uint inside=0;
    uint outside=0;

    for(auto & proj : projectionResults){
        math::Point3 projPos = math::transform(proj.transformation, position);
        if(params_.method==Method::PARITY_COUNT){
            uint parity=0;
            if(projPos[0]<proj.buffer.size.width
                    && projPos[1]<proj.buffer.size.height
                    && projPos[0]>0 && projPos[1]>0){
                //std::cout<<proj.buffer.data[projPos[0]][projPos[1]].size()<<std::endl;
                for(auto dit = proj.buffer.begin(projPos[0],projPos[1]);
                        dit != proj.buffer.end(projPos[0],projPos[1]); ++dit){

                    if(*dit<projPos[2]){
                        parity++;
                    }
                    else{
                        break;
                    }
                }
                if(parity%2==1){
                    inside++;
                    continue;
                }
                outside++;
            }
        }

        if(params_.method==Method::RAY_STABING){
            if( projPos[2]<*proj.buffer.begin(projPos[0],projPos[1])
               && projPos[2]>*proj.buffer.end(projPos[0],projPos[1])){
                inside++;
                continue;
            }
            outside++;
        }
    }

    if(params_.method==Method::PARITY_COUNT){
        if(inside>outside)
            return true;
    }
    if(params_.method==Method::RAY_STABING){
        if(outside==0)
            return true;
    }
    return false;
}

void MeshVoxelizer::visualizeDepthMap( const ProjectionResult &proj
                                      , const math::Extents3 & extents
                                      , const fs::path & path){
    cv::Mat depthMapImg( proj.buffer.size.height, proj.buffer.size.width
                     , CV_8UC3, cv::Scalar(0,0,0));

    float min = FLT_MAX;
    float max = -FLT_MAX;

    for(int x = 0; x < proj.buffer.size.width; ++x){
        for(int y = 0; y < proj.buffer.size.height; ++y){
            if(proj.buffer.begin(x,y) == proj.buffer.end(x,y)){
                continue;
            }
            min = std::min(*proj.buffer.begin(x,y),min);
            max = std::max(*std::prev(proj.buffer.end(x,y)),max);
        }
    }

    float range = max-min;

    LOG( debug )<<"Minimal/maximal depth: "<< min <<" / "<< max
                << " Range: "<<range;

    for(int x = 0; x < proj.buffer.size.width; ++x){
        for(int y = 0; y < proj.buffer.size.height; ++y){
            if(proj.buffer.begin(x,y) == proj.buffer.end(x,y)){
                depthMapImg.at<cv::Vec3b>(y,x)=
                        cv::Vec3b(0,0,255);
                continue;
            }
            float val  = (*proj.buffer.begin(x,y) - min) / range;
            depthMapImg.at<cv::Vec3b>(y,x)=
                    cv::Vec3b(255-val*255,255-val*255,255-val*255);
        }
    }

    math::Points3 extPoints;
    extPoints.push_back(math::Point3(extents.ll[0],extents.ll[1],extents.ll[2]));
    extPoints.push_back(math::Point3(extents.ll[0],extents.ll[1],extents.ur[2]));
    extPoints.push_back(math::Point3(extents.ll[0],extents.ur[1],extents.ll[2]));
    extPoints.push_back(math::Point3(extents.ll[0],extents.ur[1],extents.ur[2]));
    extPoints.push_back(math::Point3(extents.ur[0],extents.ll[1],extents.ll[2]));
    extPoints.push_back(math::Point3(extents.ur[0],extents.ll[1],extents.ur[2]));
    extPoints.push_back(math::Point3(extents.ur[0],extents.ur[1],extents.ll[2]));
    extPoints.push_back(math::Point3(extents.ur[0],extents.ur[1],extents.ur[2]));

    for(auto &point : extPoints){
        point = transform(proj.transformation,point);
    }

    cv::line(depthMapImg, cv::Point2f(extPoints[0][0],extPoints[0][1])
            , cv::Point2f(extPoints[1][0],extPoints[1][1]),cv::Scalar(255,0,0),2);
    cv::line(depthMapImg, cv::Point2f(extPoints[2][0],extPoints[2][1])
            , cv::Point2f(extPoints[3][0],extPoints[3][1]),cv::Scalar(255,0,0),2);

    cv::line(depthMapImg, cv::Point2f(extPoints[0][0],extPoints[0][1])
            , cv::Point2f(extPoints[2][0],extPoints[2][1]),cv::Scalar(255,0,0),2);
    cv::line(depthMapImg, cv::Point2f(extPoints[1][0],extPoints[1][1])
            , cv::Point2f(extPoints[3][0],extPoints[3][1]),cv::Scalar(255,0,0),2);

    cv::line(depthMapImg, cv::Point2f(extPoints[4][0],extPoints[4][1])
            , cv::Point2f(extPoints[5][0],extPoints[5][1]),cv::Scalar(255,0,0),2);
    cv::line(depthMapImg, cv::Point2f(extPoints[6][0],extPoints[6][1])
            , cv::Point2f(extPoints[7][0],extPoints[7][1]),cv::Scalar(255,0,0),2);

    cv::line(depthMapImg, cv::Point2f(extPoints[4][0],extPoints[4][1])
            , cv::Point2f(extPoints[6][0],extPoints[6][1]),cv::Scalar(255,0,0),2);
    cv::line(depthMapImg, cv::Point2f(extPoints[5][0],extPoints[5][1])
            , cv::Point2f(extPoints[7][0],extPoints[7][1]),cv::Scalar(255,0,0),2);

    cv::line(depthMapImg, cv::Point2f(extPoints[1][0],extPoints[1][1])
            , cv::Point2f(extPoints[5][0],extPoints[5][1]),cv::Scalar(255,0,0),2);
    cv::line(depthMapImg, cv::Point2f(extPoints[0][0],extPoints[0][1])
            , cv::Point2f(extPoints[4][0],extPoints[4][1]),cv::Scalar(255,0,0),2);

    cv::line(depthMapImg, cv::Point2f(extPoints[3][0],extPoints[3][1])
            , cv::Point2f(extPoints[7][0],extPoints[7][1]),cv::Scalar(255,0,0),2);
    cv::line(depthMapImg, cv::Point2f(extPoints[2][0],extPoints[2][1])
            , cv::Point2f(extPoints[6][0],extPoints[6][1]),cv::Scalar(255,0,0),2);

    cv::imwrite( path.string().c_str(), depthMapImg );
}

} //namespace geometry

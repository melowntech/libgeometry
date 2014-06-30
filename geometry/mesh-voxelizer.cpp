#ifdef _OPENMP
# include <omp.h>
#else
# define omp_get_max_threads() 1
# define omp_get_thread_num() 0
#endif

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
    if(params_.addFloor){
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

    //volume size
    math::Point3 extentsSize = extents.ur-extents.ll;
    math::Point3i volumeGridRes(
                std::ceil(extentsSize(0)/params_.voxelSize)
                ,std::ceil(extentsSize(1)/params_.voxelSize)
                ,std::ceil(extentsSize(2)/params_.voxelSize));
    math::Point3 newExtentsSize = volumeGridRes*params_.voxelSize;
    math::Point3 extentsCenter = math::center(extents);

    extents = math::Extents3( extentsCenter-newExtentsSize/2
                             , extentsCenter+newExtentsSize/2);
    LOG( info2 )<<"Volumetric grid resolution: "<<volumeGridRes;
    LOG( info2 )<<"Volume extents: "<<extents;

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
    for( auto &proj : projections){
        fs::path path = fs::path(std::string("depth")
                                 +boost::lexical_cast<std::string>(i++)
                                 +std::string(".png"));
        visualizeDepthMap(proj, extents ,path);
    }
#endif



    //create new volume
    volume_ = std::unique_ptr<Volume>(
        new Volume( extents.ll, extents.ur, params_.voxelSize, 0));

    math::Size3i vSize = volume_->cSize();
    long volMem = (long)vSize.width * vSize.height
                    * vSize.depth * sizeof(unsigned short);
    LOG( info2 )<<"Memory consumption of volume: "
        << volMem/1024.0/1024.0/1024.0 << " GB.";

    long mem = 0;
    math::Point2 avgSize(0,0);
    for(auto res : results){
        mem += res.buffer.mem();
        avgSize += math::Point2(res.buffer.size.width, res.buffer.size.height);
    }
    avgSize=avgSize/results.size();
    LOG( info2 )<<"Avg zbuffer size: "<<avgSize;
    LOG( info2 )<<"Mem zbuffer constants: "
        <<avgSize(0)/volumeGridRes(0)<<":"<<avgSize(1)/volumeGridRes(1);

    LOG( info2 )<<"Memory consumption of zbuffers: "
        << mem/1024.0/1024.0/1024.0 << " GB.";

    LOG( info2 )<<"Memory consumption of zbuffers per (m^2/gsdStr^2) "
        << mem/(volumeGridRes(0)*volumeGridRes(1))/1024.0 << " kB.";

    LOG( info2 )<<"Memory consumption of volume per (m^2/gsdStr^2)"
        << volMem/(volumeGridRes(0)*volumeGridRes(1))/1024.0 << " kB.";

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
                        volume_->set( x, y, z
                                     ,std::numeric_limits<unsigned short>::max());
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

    LOG( info2 )<<"Filling volume from seal";
    fillVolumeFromSeal();
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
    std::size_t indicesStart = 0;

    seal.vertices.push_back(
                math::Point3(extents.ll(0)-offset,extents.ll(1)-offset
                             ,extents.ll(2)-offset));
    seal.vertices.push_back(
                math::Point3(extents.ll(0)-offset,extents.ur(1)+offset
                             ,extents.ll(2)-offset));
    seal.vertices.push_back(
                math::Point3(extents.ur(0)+offset,extents.ll(1)-offset
                             ,extents.ll(2)-offset));
    seal.vertices.push_back(
                math::Point3(extents.ur(0)+offset,extents.ur(1)+offset
                             ,extents.ll(2)-offset));

    seal.addFace(indicesStart,indicesStart+2,indicesStart+1);
    seal.addFace(indicesStart+1,indicesStart+2,indicesStart+3);

    return seal;
}

void MeshVoxelizer::fillVolumeFromSeal(){
    math::Size3i vSize = volume_->cSize();
    for(int x = 0;x<vSize.width; ++x){
        for(int y = 0;y<vSize.height; ++y){
            //find first full voxel from the direction of seal
            int full = -1;
            for(int z = 0;z<vSize.depth; ++z){
                if(volume_->get(x,y,z)>0){
                        full = z;
                        break;
                }
            }
            //fill from bottom up
            for(int z = 0;z<full; ++z){
                volume_->set(x,y,z,1);
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

    uint projId = 0;
    for(auto & proj : projectionResults){
        math::Point3 projPos = math::transform(proj.transformation, position);
        projId++;
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
/*
void MeshVoxelizer::visualizeDepthMap( const ProjectionResult &proj
                                      , const math::Extents3 & extents
                                      , const fs::path & path){
    cv::Mat depthMapImg( proj.buffer.size.height, proj.buffer.size.width
                     , CV_8UC3, cv::Scalar(0,0,0));

    float min = FLT_MAX;
    float max = -FLT_MAX;

    for(int x = 0; x < proj.buffer.size.width; ++x){
        for(int y = 0; y < proj.buffer.size.height; ++y){
            if(proj.buffer.data[x][y].size()==0){
                continue;
            }
            min = std::min(proj.buffer.data[x][y].front(),min);
            max = std::max(proj.buffer.data[x][y].back(),max);
        }
    }

    float range = max-min;

    LOG( debug )<<"Minimal/maximal depth: "<< min <<" / "<< max
                << " Range: "<<range;

    for(int x = 0; x < proj.buffer.size.width; ++x){
        for(int y = 0; y < proj.buffer.size.height; ++y){
            if(proj.buffer.data[x][y].size()==0){
                depthMapImg.at<cv::Vec3b>(y,x)=
                        cv::Vec3b(0,0,255);
                continue;
            }
            float val  = (proj.buffer.data[x][y].front() - min) / range;
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

    cv::imwrite( path.native().c_str(), depthMapImg );
}
*/
} //namespace geometry

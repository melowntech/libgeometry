#include "geometry/mesh-voxelizer.hpp"
#include "math/geometry_core.hpp"
#include "imgproc/scanconversion.hpp"

#include <opencv2/opencv.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/vector.hpp>


#define RASTERIZE_MESH_DEBUG



namespace geometry{

MeshVoxelizer::MeshVoxelizer(const Parameters &params){
    params_ = params;
}

void MeshVoxelizer::add( geometry::Mesh & mesh ){
    meshes.push_back(mesh);
}

void MeshVoxelizer::voxelize(){
    if(meshes.size()<=0){
        LOG( warn2 )<<"Zero meshes to voxelize. Skipping voxelization.";
        return;
    }

    //add floor to mimic closed mesh
    if(params_.addFloor){
        for(auto & mesh: meshes){
            addSealToMesh(mesh);
        }
    }

    math::Extents3 extents(math::InvalidExtents{});

    for(auto & mesh: meshes){
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
    volume_ = std::unique_ptr<ScalarField_t<float>>(
        new ScalarField_t<float>( extents.ll, extents.ur, params_.voxelSize, 0));

    //perform rasterization in several direction
    std::vector<math::Point3> directions;
    //directions along three main axes
    directions.push_back(-math::Point3(1,0,0));
    directions.push_back(-math::Point3(0,1,0));
    directions.push_back(-math::Point3(0,0,1));

    /*
    //directions of icosahedron
    float gr = (1 + std::sqrt(5))/2;
    directions.push_back(-math::Point3((1+gr)/3,(1+gr)/3,(1+gr)/3));
    directions.push_back(-math::Point3((1+gr)/3,-(1+gr)/3,(1+gr)/3));
    directions.push_back(-math::Point3(-(1+gr)/3,(1+gr)/3,(1+gr)/3));
    directions.push_back(-math::Point3(-(1+gr)/3,-(1+gr)/3,(1+gr)/3));
    directions.push_back(-math::Point3((2+gr)/3,0,1/3.0));
    directions.push_back(-math::Point3(-(2+gr)/3,0,1/3.0));
    directions.push_back(-math::Point3(0,1/3.0,(2+gr)/3));
    directions.push_back(-math::Point3(0,-1/3.0,(2+gr)/3));
    directions.push_back(-math::Point3(1/3.0,(2+gr)/3,0));
    directions.push_back(-math::Point3(1/3.0,-(2+gr)/3,0));
    */

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


    for(const auto & dir: directions){
        math::Size2 size;
        math::Matrix4 transformation = orthoProjMat( dir,extents
                                                    , params_.voxelSize,size);
        ProjectionResult pr(transformation,LayeredZBuffer(size));
        projections.push_back(pr);
    }

    LOG( info2 )<<"Rasterizing meshes";

    //rasterize mesh using generated projection matrices
    for(auto & proj: projections){
        for(auto & mesh: meshes){
            rasterizeMesh(mesh,proj.transformation, proj.buffer);
        }
    }

    LOG( info2 )<<"Sorting buffers";
    //sort all the layered z buffers
    for( auto &proj : projections){
        for( auto &col : proj.buffer.data){
            for( auto &cell : col){
                cell.sort();
            }
        }
    }
    LOG( info2 )<<"Buffers sorted";

#ifdef RASTERIZE_MESH_DEBUG
    uint i=0;
    for( auto &proj : projections){
        fs::path path = fs::path(std::string("depth")
                                 +boost::lexical_cast<std::string>(i++)
                                 +std::string(".png"));
        visualizeDepthMap(proj, extents ,path);
    }
#endif
    uint progress = 0;
    LOG( info2 )<<"Voxelization progress: "<<progress;
    for(int x=0; x< volume_->sizeX(); ++x){
        for(int y=0; y< volume_->sizeY(); ++y){
            for(int z=0; z< volume_->sizeZ(); ++z){
                //voxel center position
                auto voxCenter
                        = volume_->grid2geo(math::Point3(x,y,z));

                if(isInside({voxCenter.x, voxCenter.y, voxCenter.z}, projections)){
                    volume_->set(x,y,z,1);
                }
            }
        }
        //voxelization progress
        uint newProgress = std::ceil((float)x/volume_->sizeX() * 100);
        if(newProgress>progress){
            progress=newProgress;
            LOG( info2 )<<"Voxelization progress: "<<progress;
        }
    }

    LOG( info2 )<<"Filling volume from seal";
    fillVolumeFromSeal();
}

std::shared_ptr<ScalarField_t<float>> MeshVoxelizer::volume(){
    return volume_;
}

void MeshVoxelizer::reset(){
    projections.clear();
}


math::Matrix4 MeshVoxelizer::orthoProjMat( const math::Point3 &direction
                            , const math::Extents3 &extents
                            , const float &voxelSize
                            , math::Size2 &viewport){
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

    LOG( debug ) << "Up" << up;
    LOG( debug ) << "Right" << right;
    LOG( debug ) << "Forward" << normalizedDir;

    ublas::subrange(c1,0,3) = right;
    ublas::subrange(c2,0,3) = up;
    ublas::subrange(c3,0,3) = normalizedDir;

    projMat = prod(rotation,projMat);

    LOG( debug )<<"Projection matrix - after rotation: "<<projMat;

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
        LOG( debug ) << point <<" - "<< pointImg;
    }

    //scale scene so the extents fit to space < -1,1 > CC/NDC
    math::Matrix4 scaleMat = ublas::identity_matrix<double>(4);
    scaleMat(0,0) = 2./size(viewExtents).width;
    scaleMat(1,1) = 2./size(viewExtents).height;

    projMat = prod(scaleMat, projMat);

    LOG( debug )<<"Projection matrix - after CC scaling: "<<projMat;

    //scale and transform scene to VC
    math::Matrix4 transformMat = ublas::identity_matrix<double>(4);
    transformMat(0,3)=1;
    transformMat(1,3)=1;

    //calculate viewport size
    LOG( debug ) << "View extents: "<< viewExtents;

    viewport.width = std::ceil(size(viewExtents).width/voxelSize);
    viewport.height = std::ceil(size(viewExtents).height/voxelSize);

    LOG( debug ) << "Viewport:" << viewport;

    scaleMat = ublas::identity_matrix<double>(4);
    scaleMat(0,0) = std::ceil(viewport.width/2);
    scaleMat(1,1) = std::ceil(viewport.height/2);

    projMat = prod(transformMat, projMat);
    projMat = prod(scaleMat, projMat);

    LOG( debug )<<"Projection matrix - final from WC to CC: "<<projMat;


    for(const auto &point : extPoints){
        math::Point3 pointImg = transform(projMat,point);
        LOG( debug ) << point <<" - "<< pointImg;
    }

    return projMat;
}

void MeshVoxelizer::addSealToMesh(geometry::Mesh & mesh){

    math::Extents3 extents = math::computeExtents(mesh.vertices);

    float offset = params_.voxelSize * params_.sealFactor;

    std::size_t indicesStart = mesh.vertices.size();
    mesh.vertices.push_back(
                math::Point3(extents.ll(0)-offset,extents.ll(1)-offset
                             ,extents.ll(2)-offset));
    mesh.vertices.push_back(
                math::Point3(extents.ll(0)-offset,extents.ur(1)+offset
                             ,extents.ll(2)-offset));
    mesh.vertices.push_back(
                math::Point3(extents.ur(0)+offset,extents.ll(1)-offset
                             ,extents.ll(2)-offset));
    mesh.vertices.push_back(
                math::Point3(extents.ur(0)+offset,extents.ur(1)+offset
                             ,extents.ll(2)-offset));

    mesh.addFace(indicesStart,indicesStart+2,indicesStart+1);
    mesh.addFace(indicesStart+1,indicesStart+2,indicesStart+3);
}

void MeshVoxelizer::fillVolumeFromSeal(){

    for(int x = 0;x<volume_->sizeX(); ++x){
        for(int y = 0;y<volume_->sizeY(); ++y){
            //find first full voxel from the direction of seal
            int full = -1;
            for(int z = 0;z<volume_->sizeZ(); ++z){
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
                      , const ProjectionResults & projectionResults){
    uint inside=0;
    uint outside=0;

    for(const auto & proj : projectionResults){
        math::Point3 projPos = math::transform(proj.transformation, position);

        if(params_.method==Method::PARITY_COUNT){
            uint parity=0;
            if(projPos[0]<proj.buffer.data.size()
                    && projPos[1]<proj.buffer.data[projPos[0]].size())
            for(const auto depth : proj.buffer.data[projPos[0]][projPos[1]]){
                if(depth<projPos[2]){
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
        if(params_.method==Method::RAY_STABING){
            if( projPos[2]<proj.buffer.data[projPos[0]][projPos[1]].back()
               && projPos[2]>proj.buffer.data[projPos[0]][projPos[1]].front()){
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

} //namespace geometry

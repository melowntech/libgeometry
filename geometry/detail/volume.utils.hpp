/**
 * @file detail/volume.utils.hpp
 * @author Tomas Drinovsky <tomas.drinovsky@citationtech.net>
 *
 * Mostly debug utilities for geometry::volume
 * Do not include directly!
 */

#ifndef GEOMETRY_VOLUME_HPP_UTILS_
#    error Do not include this implementation header directly!
#endif

/* additional functions */
template<typename Value_t>
void saveSliceAsImg( ScalarField_t<Value_t> &volume
                     , const boost::filesystem::path &path
                     , SliceDirection dir, int slice
                     , float min, float max){
    cv::Mat sliceImg;
    VolumeBase_t::Position_s pos;
    VolumeBase_t::Displacement_s dis1;
    VolumeBase_t::Displacement_s dis2;

    float range = max-min;

    switch(dir){
    case SliceDirection::X:
        sliceImg = cv::Mat( volume.sizeZ(), volume.sizeY()
                     , CV_8UC3, cv::Scalar(0,0,0));
        pos = VolumeBase_t::Position_s(slice,0,0);
        dis1 = VolumeBase_t::Displacement_s(0,1,0);
        dis2 = VolumeBase_t::Displacement_s(0,0,1);
        break;
    case SliceDirection::Y:
        sliceImg = cv::Mat( volume.sizeZ(), volume.sizeX()
                     , CV_8UC3, cv::Scalar(0,0,0));
        pos = VolumeBase_t::Position_s(0,slice,0);
        dis1 = VolumeBase_t::Displacement_s(1,0,0);
        dis2 = VolumeBase_t::Displacement_s(0,0,1);
        break;
    case SliceDirection::Z:
        sliceImg = cv::Mat( volume.sizeY(), volume.sizeX()
                     , CV_8UC3, cv::Scalar(0,0,0));
        pos = VolumeBase_t::Position_s(0,0,slice);
        dis1 = VolumeBase_t::Displacement_s(1,0,0);
        dis2 = VolumeBase_t::Displacement_s(0,1,0);
        break;
    }

    typename ScalarField_t<Value_t>::Giterator_t sit1( volume, pos, dis1 );
    typename ScalarField_t<Value_t>::Giterator_t send1 = volume.gend( sit1 );

    uint colSize = send1-sit1;

    for ( uint c = 0; c < colSize; c++ ) {
        typename ScalarField_t<Value_t>::Giterator_t sit2( volume, sit1.pos, dis2 );
        typename ScalarField_t<Value_t>::Giterator_t send2 = volume.gend( sit2 );

        uint rowSize = send2-sit2;
        for ( uint r = 0; r < rowSize; r++ ) {

            float val = std::min(std::max(sit2.value(),min),max);
            val = (val-min)/range;
            sliceImg.at<cv::Vec3b>(rowSize-1-r,c)
                    = cv::Vec3b(val*255,val*255,val*255);
            ++sit2;
        }
        ++sit1;
    }

    cv::imwrite( path.native().c_str(), sliceImg );
}

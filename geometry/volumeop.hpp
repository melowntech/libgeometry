/**
 * @file volumeop.hpp
 * @author Tomas Drinovsky <tomas.drinovsky@citationtech.net>
 *
 * Set of function and classes for various volumetric operations, including
 * filtering, morphologic operations etc.
 *
 */

#ifndef geometry_volumeop_hpp_included_
#define geometry_volumeop_hpp_included_

#include "volume.hpp"
#include "math/geometry_core.hpp"
#include "utility/enum-io.hpp"

namespace geometry{

//morphological operations

class StructuralElement{
public:
    enum class SEShape{ CUBE, SPHERE };

    StructuralElement(SEShape shape, float size);

    const std::vector<math::Point3i> & offsets() const { return offsets_; }

private:
    std::vector<math::Point3i> offsets_;
};


template <class Value_t, class Container_t>
void erosion( Container_t & container
    , const StructuralElement & se
    , const VolumeUnit<Value_t> volumeUnit = VolumeUnit<Value_t>()){

    (void) volumeUnit;
    Container_t result( container.sizeX(), container.sizeY(), container.sizeZ()
                      , volumeUnit.empty());

    utility::Progress progress(container.sizeZ());
    //for each cell in the volume
#ifdef _OPENMP
    #pragma omp parallel for schedule( dynamic, 5 )
#endif
    for(int z=0; z< container.sizeZ(); ++z){
        for(int y=0; y< container.sizeY(); ++y){
            for(int x=0; x< container.sizeX(); ++x){
                //find out if in the radius is empty voxel
                bool set = false;
                //only full pixels can be eroded
                if(container.get(x,y,z)>=volumeUnit.middle()){
                    for(const auto & offset: se.offsets()){
                        //take into account only values inside volume
                        math::Point3i lookup(x+offset(0),y+offset(1),z+offset(2));
                        if( lookup(0) >= 0 && lookup(0) < container.sizeX()
                            && lookup(1) >= 0 && lookup(1) < container.sizeY()
                            && lookup(2) >= 0 && lookup(2) < container.sizeZ()
                            && container.get(lookup(0),lookup(1),lookup(2))
                            < volumeUnit.middle()){

                            result.set(x,y,z,volumeUnit.empty());
                            set = true;
                            break;
                        }
                    }
                    if(!set)
                        result.set(x,y,z,volumeUnit.full());
                }
            }
        }
        progress.incrementAndReport(0.01);
    }
    container=std::move(result);
}

template <class Value_t, class Container_t>
void dilatation( Container_t & container
    , const StructuralElement & se
    , const VolumeUnit<Value_t> volumeUnit = VolumeUnit<Value_t>()){

    (void) volumeUnit;
    Container_t result( container.sizeX(), container.sizeY(), container.sizeZ()
                      , volumeUnit.empty());

    utility::Progress progress(container.sizeZ());
    //for each cell in the volume
#ifdef _OPENMP
    #pragma omp parallel for schedule( dynamic, 5 )
#endif
    for(int z=0; z< container.sizeZ(); ++z){
        for(int y=0; y< container.sizeY(); ++y){
            for(int x=0; x< container.sizeX(); ++x){
                //find out if in the radius is empty voxel
                bool set = false;
                //only empty pixels can become full
                if(container.get(x,y,z)<volumeUnit.middle()){
                    for(const auto & offset: se.offsets()){
                        //take into account only values inside volume
                        math::Point3i lookup(x+offset(0),y+offset(1),z+offset(2));
                        if( lookup(0) >= 0 && lookup(0) < container.sizeX()
                            && lookup(1) >= 0 && lookup(1) < container.sizeY()
                            && lookup(2) >= 0 && lookup(2) < container.sizeZ()
                            && container.get(lookup(0),lookup(1),lookup(2))
                            >= volumeUnit.middle()){

                            result.set(x,y,z,volumeUnit.full());
                            set = true;
                            break;
                        }
                    }
                    if(!set)
                        result.set(x,y,z,volumeUnit.empty());
                }
                else{
                    result.set(x,y,z,volumeUnit.full());
                }
            }
        }
        progress.incrementAndReport(0.01);
    }
    container=std::move(result);
}

UTILITY_GENERATE_ENUM_IO(StructuralElement::SEShape,
                         ((CUBE)("cube"))
                         ((SPHERE)("sphere"))
                         )



//filtering
template <class Container_t>
void filter(
    const math::FIRFilter_t & filter,
    const VolumeBase_t::Displacement_s & diff,
    Container_t & container,
    Container_t & dstVolume ) {
    typedef typename Container_t::ValueType ValueType;
    typedef Giterator_t<ValueType,Container_t> Giterator;

    assert( container.sizeX() == dstVolume.sizeX() );
    assert( container.sizeY() == dstVolume.sizeY() );
    assert( container.sizeZ() == dstVolume.sizeZ() );
    assert( diff != VolumeBase_t::Displacement_s( 0, 0, 0 ) );

    std::vector<VolumeBase_t::Position_s> poss
        = Giterator::iteratorPositions( container, diff );

    BOOST_FOREACH( VolumeBase_t::Position_s pos, poss ) {

        Giterator sit( container, pos, diff );
        Giterator send = Giterator::gend( sit );
        Giterator dit( dstVolume, pos, diff );

        int rowSize = send - sit;

        for ( int x = 0; x < rowSize; x++ ) {
            dit.setValue( filter.convolute( sit, x, rowSize ) );
            ++sit; ++dit;
        }
    }
}


template <class Container_t>
void filterInplace(
        const math::FIRFilter_t & filter,
        const VolumeBase_t::Displacement_s & diff,
        Container_t & container){

    typedef typename Container_t::ValueType ValueType;
    typedef Giterator_t<ValueType,Container_t> Giterator;

    double max=std::numeric_limits<ValueType>::max();
    double min=std::numeric_limits<ValueType>::lowest();

    std::vector<VolumeBase_t::Position_s> poss
        = Giterator::iteratorPositions( container, diff );

    for( auto pos: poss ) {

        Giterator sit( container, pos, diff );
        Giterator send = Giterator::gend( sit );

        int rowSize = send - sit;

        std::vector<ValueType> filtered(rowSize);
        auto dit = filtered.begin();

        for ( int x = 0; x < rowSize; x++ ) {
            *dit=(ValueType)std::min(
                    std::max(filter.convolute( sit, x, rowSize ),min),max);
            ++sit; ++dit;
        }

        //write filtered values into source volume
        dit = filtered.begin();
        sit.pos = pos;
        for ( int x = 0; x < rowSize; x++ )
        {
                sit.setValue(*dit);
                ++sit; ++dit;
        }
    }
}

/*
 * Specialized version of filtering used for VolumeArrays. The container
 * representation allows to use multithreding.
 */
template <class Value_t>
void filterInplace(
        const math::FIRFilter_t & filter,
        const VolumeBase_t::Displacement_s & diff,
        VolumeArray<Value_t> & container){

    typedef Value_t ValueType;
    typedef Giterator_t<Value_t,VolumeArray<Value_t>> Giterator;

    double max=std::numeric_limits<ValueType>::max();
    double min=std::numeric_limits<ValueType>::lowest();

    std::vector<VolumeBase_t::Position_s> poss
        = Giterator::iteratorPositions( container, diff );
#ifdef _OPENMP
    #pragma omp parallel for schedule( dynamic, 20 )
#endif
    for(uint p=0; p<poss.size(); ++p){
        VolumeBase_t::Position_s pos = poss[p];
        Giterator sit( container, pos, diff );
        Giterator send = Giterator::gend( sit );

        int rowSize = send - sit;

        std::vector<ValueType> filtered(rowSize);
        auto dit = filtered.begin();

        for ( int x = 0; x < rowSize; x++ ) {
            *dit=(ValueType)std::min(
                    std::max(filter.convolute( sit, x, rowSize ),min),max);
            ++sit; ++dit;
        }

        //write filtered values into source volume
        dit = filtered.begin();
        sit.pos = pos;
        for ( int x = 0; x < rowSize; x++ )
        {
                sit.setValue(*dit);
                ++sit; ++dit;
        }
    }
}

} //namespace geometry

#endif

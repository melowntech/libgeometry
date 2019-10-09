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
 * @file volumeop.hpp
 * @author Tomas Drinovsky <tomas.drinovsky@citationtech.net>
 *
 * Set of function and classes for various volumetric operations, including
 * filtering, morphologic operations etc.
 *
 */

#ifndef geometry_volumeop_hpp_included_
#define geometry_volumeop_hpp_included_

#include "geometry/volume.hpp"
#include "math/geometry_core.hpp"
#include "utility/enum-io.hpp"

namespace geometry {

//morphological operations

class StructuralElement{
public:
    enum class SEShape{ CUBE, SPHERE };

    StructuralElement(SEShape shape, float size);

    const std::vector<math::Point3i> & offsets() const { return offsets_; }

private:
    std::vector<math::Point3i> offsets_;
};


template <typename Value_t, class Container_t>
void erosion( Container_t & container
    , const StructuralElement & se
    , const VolumeUnit<Value_t> volumeUnit = VolumeUnit<Value_t>()){

    (void) volumeUnit;
    Container_t result( container.sizeX(), container.sizeY(), container.sizeZ()
                      , volumeUnit.empty());
    //for each cell in the volume
#ifdef _OPENMP
    #pragma omp parallel for schedule( dynamic, 5 )
#endif
    for(int z=0; z< container.sizeZ(); ++z){
        for(int y=0; y< container.sizeY(); ++y){
            for(int x=0; x< container.sizeX(); ++x){
                //find out if in the radius is empty voxel
                //only full pixels can be eroded
                Value_t cvalue = container.get(x,y,z);
                for(const auto & offset: se.offsets()){
                    //take into account only values inside volume
                    math::Point3i lookup(x+offset(0),y+offset(1),z+offset(2));
                    if( lookup(0) >= 0 && lookup(0) < container.sizeX()
                        && lookup(1) >= 0 && lookup(1) < container.sizeY()
                        && lookup(2) >= 0 && lookup(2) < container.sizeZ()){
                        cvalue = std::min(container.get(lookup(0),lookup(1),lookup(2)),cvalue);
                    }
                }
                result.set(x,y,z,cvalue);
            }
        }
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

    //for each cell in the volume
#ifdef _OPENMP
    #pragma omp parallel for schedule( dynamic, 5 )
#endif
    for(int z=0; z< container.sizeZ(); ++z){
        for(int y=0; y< container.sizeY(); ++y){
            for(int x=0; x< container.sizeX(); ++x){
                //find out if in the radius is empty voxel
                //only empty pixels can become full
                Value_t cvalue = container.get(x,y,z);
                for(const auto & offset: se.offsets()){
                    //take into account only values inside volume
                    math::Point3i lookup(x+offset(0),y+offset(1),z+offset(2));
                    if( lookup(0) >= 0 && lookup(0) < container.sizeX()
                        && lookup(1) >= 0 && lookup(1) < container.sizeY()
                        && lookup(2) >= 0 && lookup(2) < container.sizeZ()){
                        cvalue = std::max(container.get(lookup(0),lookup(1),lookup(2)),cvalue);
                    }
                }
                result.set(x,y,z,cvalue);
            }
        }
    }
    container=std::move(result);
}

UTILITY_GENERATE_ENUM_IO(StructuralElement::SEShape,
                         ((CUBE)("cube"))
                         ((SPHERE)("sphere"))
                         )



//filtering
/** Filters volumetric data in particular direction
 * \param filter FIR filter used for filtering
 * \param diff defines displacement used as step
 * \param container source volumetric data
 * \param dstVolume destination volumetric data
 */
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

/** Filters volumetric data in x,y,z axis
 * \param filter FIR filter used for filtering
 * \param container source volumetric data
 * \param dstVolume destination volumetric data
 */
template <class Container_t>
void filter(
    const math::FIRFilter_t & filter,
    Container_t & container,
    Container_t & dstVolume ) {

    typename VolumeBase_t::Displacement_s directions[3];
    directions[0] = typename VolumeBase_t::Displacement_s(1,0,0);
    directions[1] = typename VolumeBase_t::Displacement_s(0,1,0);
    directions[2] = typename VolumeBase_t::Displacement_s(0,0,1);

    for(uint fAxis = 0; fAxis<3; ++fAxis){
        LOG( info2 )<<"Filtering volume in axis "<<fAxis;
        filter( filter, directions[fAxis], container, dstVolume);
    }
}

/** Filters volumetric data inplace in particular direction
 * \param filter FIR filter used for filtering
 * \param diff defines displacement used as step
 * \param container source volumetric data
 */
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

    UTILITY_OMP(parallel for schedule( dynamic, 20 ))
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

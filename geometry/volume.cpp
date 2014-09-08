/*
 * Volume.cpp
 */

#include "volume.hpp"

#include <math/filters.hpp>
#include <boost/foreach.hpp>
#include "volumeop.hpp"

/**
 * Following code was untested after the volume related classes were
 * remade to use templated container. In the time of this change no part
 * of the this code were used in any part of the vadstena toolkit.
 */

namespace geometry {
/** class BitfieldReconstruction_t */

BitfieldReconstruction_t::BitfieldReconstruction_t( const Bitfield_t & from,
    const double delta, double filterCutoffPeriod  )
    : ScalarField_t<float,VolumeOctree<float>>( from.lower(), from.upper(),
                            from.voxelSize(), -1.0 ) {

    typedef Giterator_t<float,VolumeOctree<float>> Giterator;

    // obtain distance field
    LOG( info2 ) << "Obtaining distance map.";
    DistanceMap_t<float> distanceMap( from, delta / 2.0 * 1.1 );

    // initialize voting field
    LOG( info2 ) << "Creating voting field.";
    VotingField_t vfield( from );

    // set the 13 scanning direcitons
    std::vector<VolumeBase_t::Displacement_s> dspls;

    dspls.push_back( VolumeBase_t::Displacement_s( 1, 0, 0 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 0, 1, 0 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 0, 0, 1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 1, 0 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, -1, 0 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 0, 1, 1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 0, 1, -1  ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 0,  1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 0, -1  ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 1, 1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 1, -1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, -1, 1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( -1, 1, 1 ) );

    // process all scanlines
    BOOST_FOREACH( VolumeBase_t::Displacement_s diff, dspls ) {

        LOG ( info1 ) << "Processing direction " << diff;

        std::vector<VolumeBase_t::Position_s> poss
            = Giterator::iteratorPositions(distanceMap.container(), diff );

        BOOST_FOREACH( VolumeBase_t::Position_s pos, poss ) {

            Giterator begin =  Giterator::gbegin(distanceMap.container(), pos, diff );
            Giterator end =  Giterator::gend(begin);

            scanline( begin, end, vfield, delta / 2.0 );
       }
    }

    // iterate through voting field, evaluating polls
    LOG( info2 ) << "Evaluating polls.";
    ScalarField_t<float,VolumeOctree<float>> rawfield( _lower, _upper, _voxelSize, -1.0 );

    for ( int i = 0; i < container_.sizeX(); i++ )
        for (  int j = 0; j < container_.sizeY(); j++ )
            for ( int k = 0; k < container_.sizeZ(); k++ )
                rawfield.set( i, j, k, pollResult( vfield.get( i, j, k ) ) );

    // perform filtering
    LOG( info2 ) << "Filtering output.";
    math::LowPassFilter_t lPFilter( filterCutoffPeriod * delta / _voxelSize,
        uint( filterCutoffPeriod * ceil( delta / _voxelSize ) ) );

    LOG( info1 ) << "Filtering in x.";
    filter( lPFilter, Displacement_s( 1, 0, 0 )
            ,rawfield.container(), rawfield.container());

    LOG( info1 ) << "Filtering in y.";
    filter( lPFilter, Displacement_s( 0, 1, 0 )
            ,rawfield.container(), rawfield.container());

    LOG( info1 ) << "Filtering in z.";
    filter( lPFilter, Displacement_s( 0, 0, 1 )
            ,rawfield.container(), rawfield.container());

    // all done
}

BitfieldReconstruction_t::BitfieldReconstruction_t( const PointCloud & cloud,
        double voxelSize,
        double delta,
        double filterCutoffPeriod )
    : ScalarField_t<float,VolumeOctree<float>>( cloud.lower(), cloud.upper(),
                            voxelSize, -1.0 ) {

    typedef Giterator_t<float,VolumeOctree<float>> Giterator;

    // obtain distance field
    LOG( info2 ) << "Obtaining distance map.";
    DistanceMap_t<float> distanceMap( cloud, voxelSize, delta / 2.0 * 1.1 );
    _lower = distanceMap.lower(); _upper = distanceMap.upper();

    // initialize voting field
    LOG( info2 ) << "Creating voting field.";
    VotingField_t vfield( *this );

    // set the 13 scanning direcitons
    std::vector<VolumeBase_t::Displacement_s> dspls;

    dspls.push_back( VolumeBase_t::Displacement_s( 1, 0, 0 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 0, 1, 0 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 0, 0, 1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 1, 0 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, -1, 0 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 0, 1, 1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 0, 1, -1  ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 0,  1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 0, -1  ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 1, 1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, 1, -1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( 1, -1, 1 ) );
    dspls.push_back( VolumeBase_t::Displacement_s( -1, 1, 1 ) );

    // process all scanlines
    BOOST_FOREACH( VolumeBase_t::Displacement_s diff, dspls ) {

        LOG ( info1 ) << "Processing direction " << diff;


        std::vector<VolumeBase_t::Position_s> poss
            = Giterator::iteratorPositions(distanceMap.container(), diff );

        BOOST_FOREACH( VolumeBase_t::Position_s pos, poss ) {

            Giterator begin =  Giterator::gbegin(distanceMap.container(), pos, diff );
            Giterator end =  Giterator::gend(begin);

            scanline( begin, end, vfield, delta / 2.0 );
       }
    }

    // iterate through voting field, evaluating polls
    LOG( info2 ) << "Evaluating polls.";
    ScalarField_t<float,VolumeOctree<float>>
        rawfield( _lower, _upper, _voxelSize, -1.0 );

    for ( int i = 0; i < container_.sizeX(); i++ )
        for (  int j = 0; j < container_.sizeY(); j++ )
            for ( int k = 0; k < container_.sizeZ(); k++ )
                rawfield.set( i, j, k, pollResult( vfield.get( i, j, k ) ) );

    // perform filtering
    LOG( info2 ) << "Filtering output.";
    math::LowPassFilter_t lPFilter( filterCutoffPeriod * delta / _voxelSize,
        uint( filterCutoffPeriod * ceil( delta / _voxelSize ) ) );

    LOG( info1 ) << "Filtering in x.";
    filter( lPFilter, Displacement_s( 1, 0, 0 )
            ,rawfield.container(), rawfield.container());

    LOG( info1 ) << "Filtering in y.";
    filter( lPFilter, Displacement_s( 0, 1, 0 )
            ,rawfield.container(), rawfield.container());

    LOG( info1 ) << "Filtering in z.";
    filter( lPFilter, Displacement_s( 0, 0, 1 )
            ,rawfield.container(), rawfield.container());

    // all done
}


void BitfieldReconstruction_t::scanline(
    const Giterator_t<float,VolumeOctree<float>> & begin,
    const Giterator_t<float,VolumeOctree<float>> & end,
    VotingField_t & vfield,
    const double delta ) {

    typedef enum { IN_N, OUT_N } State_t;

    typedef Giterator_t<float,VolumeOctree<float>> Giterator;

    std::vector<Giterator> boundaryPoints;
    State_t state = OUT_N;
    Giterator entryPoint( begin );
    int incount = 0;

    // first, guess where the boundary intersection points are
    boundaryPoints.push_back( begin );

    for ( Giterator it = begin; it < end; ++it ) {

        if ( state == OUT_N ) {
            // outside zone of interest
            if ( it.value() < delta ) {
                // entering zone of interest
                entryPoint = it;
                incount = 0;
                state = IN_N;
            }

            continue;
        }

        if ( state == IN_N ) {
            // inside zone of interest
            if ( it.value() > delta ) {
                // leaving zone of interest, evaluate
                incount >>= 1;
                boundaryPoints.push_back( entryPoint + incount );

                state = OUT_N;

            } else
                incount++;

            continue;
        }
    }

    boundaryPoints.push_back( end );

    // only scanlines with even number of boundary points count
    if ( math::odd( boundaryPoints.size() ) ) return;

    // each two consequent boundary points define a positive vote segment
    for ( uint i = 0; i <= boundaryPoints.size() - 2; i++ ) {

        Giterator bbg = boundaryPoints[i];
        Giterator ben = boundaryPoints[i+1];

        VotingField_t::Giterator vbg( vfield, bbg.pos, bbg.diff );
        VotingField_t::Giterator ven( vfield, ben.pos, ben.diff );

        for ( VotingField_t::Giterator it = vbg; it < ven; ++it ) {

            Poll_s poll = it.value();
            if ( math::odd( i ) )
                poll.positives++;
            else
                poll.negatives++;
            it.setValue( poll );
        }
    }

    // all done
}

float BitfieldReconstruction_t::pollResult( const Poll_s & poll ) {

    //if ( poll.positives >= 5 ) return 1.0; else return -1.0;
    //return ( poll.positives > poll.negatives );
    if ( poll.positives == 0 )
        return -1.0;
    else
        return ( (float) poll.positives - poll.negatives )
            / ( (float) poll.positives + poll.negatives );
}

} // namespace geometry

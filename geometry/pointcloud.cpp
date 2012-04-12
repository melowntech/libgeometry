/* pointcloud.cpp */

#include "pointcloud.hpp"

#include <dbglog/dbglog.hpp>
#include <boost/foreach.hpp>
#include <fstream>

namespace geometry {


void PointCloud_t::push_back( const ublas::vector<double> & x ) {

    updateExtents( x );
    std::vector<ublas::vector<double> >::push_back( x );
}

PointCloud_t::iterator PointCloud_t::insert( iterator position,
    const ublas::vector<double> & x ) {

    updateExtents( x );
    return std::vector<ublas::vector<double> >::insert( position, x );
}

void PointCloud_t::insert( iterator position, size_type n,
    const ublas::vector<double> & x ) {

    updateExtents( x );
    std::vector<ublas::vector<double> >::insert( position, n, x );
}

void PointCloud_t::clear() {
    _upper = _lower = ublas::zero_vector<double>(3);
}


void PointCloud_t::dump( const std::string & path ) {

    std::fstream f;

    f.exceptions( std::ios::badbit | std::ios::failbit );

    try {
    
        f.open( path, std::ios_base::out | std::ios_base::trunc );
    
        for ( const_iterator it = begin(); it < end(); ++it ) {
            f << (*it)(0) << "\t" << (*it)(1) << "\t" << (*it)(2) << "\n";
        }

    } catch ( std::ios_base::failure & e ) {
        
        LOG( err2 ) << "Failed to write to '" << path
                    << "', error: " << e.what();
        throw;
    }
}

void PointCloud_t::updateExtents( const ublas::vector<double> & x ) {

    if ( empty() ) {

        _lower = _upper = x; return;
    }

    if ( x[0] < _lower[0] ) _lower[0] = x[0];
    if ( x[1] < _lower[1] ) _lower[1] = x[1];
    if ( x[2] < _lower[2] ) _lower[2] = x[2];
    if ( x[0] > _upper[0] ) _upper[0] = x[0];
    if ( x[1] > _upper[1] ) _upper[1] = x[1];
    if ( x[2] > _upper[2] ) _upper[2] = x[2];
}

double PointCloud_t::samplingDelta( float bulkThreshold ) const {

    ThreeDistance_t * distArray = new ThreeDistance_t[ size() ];

    // sanity
    assert( ! empty() );

    // obtain array of closest neighbour distances
    double maxDist = ublas::norm_2( _upper - _lower );
    for ( uint i = 0; i < size(); i++ )
        distArray[i] = ThreeDistance_t( maxDist );

    for ( uint i = 0; i < size(); i++ )
        for ( uint j = 0; j < i; j++ ) {

            distArray[i].update( at( i ) - at( j ) );
            distArray[j].update( at( i ) - at( j ) );
       }

    // sort array
    std::sort( distArray, distArray + size() );

    // return delta
    double retval =  distArray[ int( floor( bulkThreshold * ( size() - 1 ) ) ) ].value();
    delete[] distArray;
    return retval;
}

/* PointCloud_t::ThreeDistance_s */

void PointCloud_t::ThreeDistance_t::update( const ublas::vector<double> diff ) {

    double dist = ublas::norm_2( diff );

    if ( dist < VERY_SMALL_NUMBER ) return;

    int code = ( fabs( diff[0] ) > fabs( diff[1] ) ) << 0
    | ( fabs( diff[1] ) > fabs( diff[2] ) ) << 1
    | ( fabs( diff[2] ) > fabs( diff[0] ) ) << 2;

    if ( code == 0 || code == 1 || code == 3 ) {

        if ( dist < distX ) distX = dist;
    }

    if ( code == 2 || code == 6 ) {

        if ( dist < distY ) distY = dist;
    }

    if ( code == 4 || code == 5 ) {

        if ( dist < distZ ) distZ = dist;
    }

}


double PointCloud_t::ThreeDistance_t::value() const {

    double dists[3];
    dists[0] = distX; dists[1] = distY; dists[2] = distZ;
    std::sort( dists, dists + 3 );
    return dists[ 1 ];
}

} // namespace geometry

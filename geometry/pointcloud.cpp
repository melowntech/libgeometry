/* pointcloud.cpp */

#include "pointcloud.hpp"

#include <dbglog/dbglog.hpp>
#include <boost/foreach.hpp>
#include <fstream>

namespace geometry {


void PointCloud::push_back( const math::Point3 & x ) {

    updateExtents( x );
    std::vector<math::Point3>::push_back( x );
}

PointCloud::iterator PointCloud::insert( iterator position,
    const math::Point3 & x ) {

    updateExtents( x );
    return std::vector<math::Point3>::insert( position, x );
}

void PointCloud::insert( iterator position, size_type n,
    const math::Point3 & x ) {

    updateExtents( x );
    std::vector<math::Point3>::insert( position, n, x );
}

void PointCloud::clear() {
    _upper = _lower = ublas::zero_vector<double>(3);
}


void PointCloud::dump( const std::string & path ) {

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

void PointCloud::updateExtents( const math::Point3 & x ) {

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

double PointCloud::samplingDelta( float bulkThreshold ) const {

    ThreeDistance * distArray = new ThreeDistance[ size() ];

    // sanity
    assert( ! empty() );

    // obtain array of closest neighbour distances
    double maxDist = ublas::norm_2( _upper - _lower );
    for ( uint i = 0; i < size(); i++ )
        distArray[i] = ThreeDistance( maxDist );

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

/* PointCloud::ThreeDistance */

void PointCloud::ThreeDistance::update( const math::Point3 diff ) {

    double dist = ublas::norm_2( diff );

    if ( dist < 1E-15 ) return;

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


double PointCloud::ThreeDistance::value() const {

    double dists[3];
    dists[0] = distX; dists[1] = distY; dists[2] = distZ;
    std::sort( dists, dists + 3 );
    return dists[ 1 ];
}

} // namespace geometry

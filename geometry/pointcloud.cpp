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
    extents_ = math::Extents3();
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


void PointCloud::load( const std::string & path ) {

    std::ifstream f;
    f.exceptions( std::ios::badbit );

    try {

        f.open( path );
        clear();

        double x, y, z;
        while (f >> x >> y >> z) {
            push_back( math::Point3(x, y, z) );
        }

    } catch ( std::ios_base::failure & e ) {

        LOG( err2 ) << "Failed to read '" << path
                    << "', error: " << e.what();
        throw;
    }
}


void PointCloud::updateExtents( const math::Point3 & x ) {
    if ( empty() ) {
        extents_ = math::Extents3(x, x);
        return;
    }

    update(extents_, x);
}

double PointCloud::samplingDelta( float bulkThreshold ) const {

    ThreeDistance * distArray = new ThreeDistance[ size() ];

    // sanity
    assert( ! empty() );

    // obtain array of closest neighbour distances
    double maxDist = ublas::norm_2( extents_.ur - extents_.ll );
    for ( unsigned int i = 0; i < size(); i++ )
        distArray[i] = ThreeDistance( maxDist );

    for ( unsigned int i = 0; i < size(); i++ )
        for ( unsigned int j = 0; j < i; j++ ) {

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

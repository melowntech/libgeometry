/**
 * @file pointcloud.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Point clouds, or set of 3D points (usually surface boundary samples).
 */

#ifndef GEOMETRY_POINTCLOUD_HPP
#define GEOMETRY_POINTCLOUD_HPP

#include <set>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

namespace geometry {

/**
 * PointCloud is essentially a std::vector of 3D points, with extents
 * maintenance and sampling density computation.
 */

class PointCloud_t: public std::vector<ublas::vector<double> > {

public :

    /** Initialize */
    PointCloud_t() :
        _lower( ublas::zero_vector<double>( 3 ) ),
        _upper( ublas::zero_vector<double>( 3 ) ) {};

    /** push_back */
    void push_back ( const ublas::vector<double> & x );

    /** insert */
    iterator insert ( iterator position, const ublas::vector<double> & x );

    /** insert */
    void insert ( iterator position, size_type n,
                  const ublas::vector<double> & x );

    /** insert */
    template <class InputIterator>
    void insert ( iterator position, InputIterator first, InputIterator last );

    /** clear */
    void clear();

    /**
     * Save to a file. The format is simplistic, with one line per point,
     * three whitespace separated values per line.
     */
    void dump( const std::string & path );

    /**
     * Return a measure of euclidian distance to a nearest point.
     * Sampling delta * 100 % point are at most as far from a nearest point as
     * the return value.
     */
    double samplingDelta ( float bulkThreshold = 0.5 ) const;

    /** Upper bound of all points */
    ublas::vector<double> upper() const { assert( ! empty() ); return _upper; }

    /** Upper bound of all points. */
    ublas::vector<double> lower() const { assert( ! empty() ); return _lower; }

private :

    /* forbidden modifiers */
    template <class InputIterator>
    void assign ( InputIterator first, InputIterator last ) {
        assert( false ); }
    void pop_back ( ) { assert( false); }
    iterator erase ( iterator position ) { (void) position; assert( false ); }
    iterator erase ( iterator first, iterator last ) {
        (void) first; (void) last;
        assert( false );
    }
    void swap( std::vector<ublas::vector< double> >& vec ) {
        (void) vec;
        assert( false ); }

    void updateExtents( const ublas::vector<double> & x );

    class ThreeDistance_t {

    public :

        ThreeDistance_t( double value = 0.0 )
            : distX( value ), distY( value ), distZ( value ) {}

        void update( const ublas::vector<double> diff );

        double value() const;

        bool operator < ( const ThreeDistance_t & op ) const {

            return value() < op.value();
        }
    private :

        double distX, distY, distZ;
    };

    ublas::vector<double> _lower, _upper;
};


/* template method implementation */

template <class InputIterator>
void PointCloud_t::insert ( iterator position, InputIterator first,
                            InputIterator last ) {

    for ( InputIterator it = first; it < last; it++ )
        updateExtents( *it );

    std::vector<ublas::vector<double > >::insert( position, first, last );
}


} // namespace geometry

#endif

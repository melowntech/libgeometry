/**
 * @file volume.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Generic volumetric modeling class.
 *
 * This module provides a class for generic volumetric modeling. Its core is
 * a class representing a generic quantity in a voxel grid, scalar or vector.
 * The voxel grid is represented as an octree.
 */


#ifndef GEOMETRY_VOLUME_HPP
#define GEOMETRY_VOLUME_HPP

#include <math/math_all.hpp>
#include <dbglog/dbglog.hpp>

#include <boost/foreach.hpp>
#include <set>
#include <vector>

namespace geometry {

class VolumeBase_t {

public:

    class Displacement_s;

    struct Position_s {

        int x, y, z;

        Position_s( int posx = 0, int posy = 0, int posz = 0 )
            : x( posx ), y( posy ), z( posz ) {}

        Displacement_s operator - ( const Position_s & pos2 ) const {
            return Displacement_s( x - pos2.x, y - pos2.y, z - pos2.z );
        }

        Position_s operator + ( const Displacement_s & diff ) const {
            return Position_s( x + diff.x, y + diff.y, z + diff.z );
        }

        Position_s operator - ( const Displacement_s & diff ) const {
            return Position_s( x - diff.x, y - diff.y, z - diff.z );
        }

        bool operator == ( const Position_s & op ) const {
            return( x == op.x && y == op.y && z == op.z );
        }

        bool operator < ( const Position_s & op2 ) const {
            if ( z < op2.z
                 || ( z == op2.z && y < op2.y )
                 || ( z == op2.z && y == op2.y && x < op2.x ) )
                return true;
            else
                return false;
        }
    };

    struct Displacement_s {
        int x, y, z;

        Displacement_s( int difx = 0, int dify = 0, int difz = 0 )
        : x( difx ), y( dify ), z( difz ) {}

        bool operator == ( const Displacement_s & op ) const {
            return( x == op.x && y == op.y && z == op.z );
        }

        bool operator != ( const Displacement_s & op ) const {
            return( x != op.x || y != op.y || z != op.z );
        }

        Displacement_s operator * ( const int f ) const {
            return Displacement_s( f * x, f * y, f * z );
        }
    };


    struct FPosition_s {
        double x, y, z;
        FPosition_s( const double x = 0.0, const double y = 0.0,
                     const double z = 0.0 ) : x( x ), y( y ), z( z ) {}

        FPosition_s( const ublas::vector<double> & op ) :
            x( op[0] ), y( op[1] ), z( op[2] ) {};
    };
};

template <typename E, typename T>
std::basic_ostream<E, T> & operator << (
    std::basic_ostream<E,T> & os, const VolumeBase_t::FPosition_s & pos ) {

    os << "[" << pos.x << "," << pos.y << "," << pos.z << "]";
    return os;
}

template <typename E, typename T>
std::basic_ostream<E, T> & operator << (
    std::basic_ostream<E,T> & os, const VolumeBase_t::Displacement_s & pos ) {

    os << "[" << pos.x << "," << pos.y << "," << pos.z << "]";
    return os;
}


template <typename Value_t>
class Volume_t : public VolumeBase_t {

public:


    /** Construct a volume and initialize it to a given value. */
    Volume_t( const int sizeX, const int sizeY, const int sizeZ,
              const Value_t & initValue );

    /** Volume destruction. */
    ~Volume_t();

    /** Value getter. */
    Value_t get( int i, int j, int k ) const;

    /** Value setter. */
    void set( int i, int j, int k, const Value_t & value );

    /**
     * Generic volume iterators, defined by starting position and displacement
     * vector. This is not a true iterator, as dereferencing is not possible -
     * use value() and setValue() instead.
     */
    struct Giterator_t {

        Volume_t<Value_t> * volume;
        Position_s pos;
        Displacement_s diff;

        Giterator_t( Volume_t<Value_t> & volume,
            const Position_s & pos = Position_s(),
            const Displacement_s & diff = Displacement_s() )
            : volume( & volume ), pos( pos ), diff( diff ) {}

        Giterator_t & operator++()  { pos = pos + diff; return *this; }

        Giterator_t operator + ( const int & count ) const {
            return Giterator_t( *volume, pos + diff * count, diff );
        }

        Giterator_t operator - ( const int & count ) const {
            return Giterator_t( *volume, pos - diff * count, diff );
        }

        int operator - ( const Giterator_t & op ) const {

            assert( op.diff == diff );
            Displacement_s bdiff = pos - op.pos;

            if ( diff.x > 0 )
                return bdiff.x / diff.x;

            if ( diff.y > 0 )
                return bdiff.y / diff.y;

            if ( diff.z > 0 )
                return bdiff.z / diff.z;

            return 0;
        }

        Value_t value() const { return volume->get( pos.x, pos.y, pos.z ); }

        Value_t operator[] ( int i ) const {
            return ( *this + i ).value();
        }

        void setValue( const Value_t & value ) const {
            volume->set( pos.x, pos.y, pos.z, value ); }

        bool operator < ( const Giterator_t & s );
        bool operator <= ( const Giterator_t & s );
    };

    /** Iterator initialization */
    Giterator_t gbegin( const Position_s & pos,
        const Displacement_s & diff ) {
            return Giterator_t( *this, pos, diff ); }

    /** Iterator end marker */
    Giterator_t gend( const Giterator_t & begin );

    /** Helper function used to provide, for a given displacement vector, a set
     * of iterator initial position such that iterators with such displacement
     * cover the entire volume. */
    std::set<Position_s> iteratorPositions (
        const Displacement_s & diff ) const;

    int sizeX() const { return _sizeX; }
    int sizeY() const { return _sizeY; }
    int sizeZ() const { return _sizeZ; }

protected:
    struct Node_s {
        typedef enum { SOLID, GRAY } Type_t;

        const static unsigned char OCT_X = 0x04;
        const static unsigned char OCT_Y = 0x02;
        const static unsigned char OCT_Z = 0x01;



        typedef enum {
            LBB = 0x00,
            LBF = OCT_Z,
            LTB = OCT_Y,
            LTF = OCT_Y | OCT_Z,
            RBB = OCT_X,
            RBF = OCT_X | OCT_Z,
            RTB = OCT_X | OCT_Y,
            RTF = OCT_X | OCT_Y | OCT_Z } Octant_t;

        Type_t type;
        Value_t value;
        Node_s * subnodes[8]; // indexed by quadrant

        Octant_t findOctant( int nodeSize, const Position_s & pos );

        Position_s toOctant( Octant_t octant, int nodeSize,
            const Position_s & pos );
        Position_s fromOctant( Octant_t octant, int nodeSize,
            const Position_s & pos );

        Value_t get( int nodeSize, const Position_s & pos );
        void set( int nodeSize, const Position_s & pos, const Value_t & value );

        Node_s( const Value_t & value ) : type( SOLID ), value( value ) {}

        ~Node_s();
    };

    Node_s * _root;
    int _rootSize;
    int _sizeX, _sizeY, _sizeZ;
    Value_t _initValue;
};

/** GeoVolume is a volume with defined floating point georeferencing. */
template <typename Value_t>
class GeoVolume_t : public Volume_t<Value_t> {

public :

    GeoVolume_t( const VolumeBase_t::FPosition_s & lower,
                 const VolumeBase_t::FPosition_s & upper,
                 const double voxelSize, const Value_t & initValue );

    VolumeBase_t::FPosition_s lower() const { return _lower; }
    VolumeBase_t::FPosition_s upper() const { return _upper; }
    double voxelSize() const { return _voxelSize; }

    Value_t fget( const double x, const double y, const double z );
    void fset( const double x, const double y, const double z,
        const Value_t & value );

protected :

    typename VolumeBase_t::Position_s geo2grid(
        const VolumeBase_t::FPosition_s & gpos );
    VolumeBase_t::FPosition_s grid2geo(
        const typename VolumeBase_t::Position_s & pos );
    VolumeBase_t::FPosition_s grid2geo( const VolumeBase_t::FPosition_s & pos );

    VolumeBase_t::FPosition_s _lower, _upper;
    double _voxelSize;
};

/** ScalarField is a geovolume with scalar values. */
template <typename Value_t>
class ScalarField_t : public GeoVolume_t<Value_t> {

public :

    typedef enum { TO_MIN, TO_MAX } SurfaceOrientation_t;

    ScalarField_t(
        const typename GeoVolume_t<Value_t>::FPosition_s & lower,
        const typename GeoVolume_t<Value_t>::FPosition_s & upper,
        const double voxelSize, const Value_t & initValue )
        : GeoVolume_t<Value_t>( lower, upper, voxelSize, initValue ) {};


    template <typename DstVolume_t>
    void filter( 
        const math::FIRFilter_t & filter,
        const VolumeBase_t::Displacement_s & diff,
        DstVolume_t & dstVolume );

    /**
     * Provide basic visualization of a scalar field isosurface as a set of
     * quads, separating voxels on differnt sides of the isosurface.
     * The output is a list of points, where each consequent quadruple defines
     * a quad.
     */
    std::vector<typename GeoVolume_t<Value_t>::FPosition_s>
        getQuads( const Value_t & threshold,
            const SurfaceOrientation_t orientation = TO_MIN );

    /**
     * Extract isosurface with a marching tetrahedrons algorithm.
     * The output is a list of points where each consequent triple defines a
     * triangle.
     */
    std::vector<typename VolumeBase_t::FPosition_s>
        isosurface( const Value_t & threshold,
            const SurfaceOrientation_t orientation = TO_MIN );

private:

    /** Used for isosurface extraction */
    void isoFromTetrahedron(
        std::vector<typename GeoVolume_t<Value_t>::FPosition_s> & retval,
        const typename GeoVolume_t<Value_t>::FPosition_s & vx0,
        const Value_t & value0,
        const typename GeoVolume_t<Value_t>::FPosition_s & vx1,
        const Value_t & value1,
        const typename GeoVolume_t<Value_t>::FPosition_s & vx2,
        const Value_t & value2,
        const typename GeoVolume_t<Value_t>::FPosition_s & vx3,
        const Value_t & value3,
        const Value_t & threshold, const SurfaceOrientation_t orientation );

    /** used for isosurface extraction */
    typename GeoVolume_t<Value_t>::FPosition_s interpolate(
            const typename GeoVolume_t<Value_t>::FPosition_s & p1,
            const Value_t & value1,
            const typename GeoVolume_t<Value_t>::FPosition_s & p2,
            const Value_t & value2,
            Value_t midval );

};

/** Bitfield is a geovolume with true/false values. */
typedef ScalarField_t<bool> Bitfield_t;

/** Distance map provides a distance map for a bitfield. Each element
  * of a distance map holds the euclidian distance from the nearest
  * non zero point. */

template <typename Value_t>
class DistanceMap_t: public ScalarField_t<Value_t> {
public:

    /**
     * Create distance map from a bitfield. InitValue corresponds to the
     * maximum distance (infty). The lower this value, the more efficient
     * the memory representation.
     */
    DistanceMap_t( const Bitfield_t & bitfield, const Value_t initValue );



private:
    struct DistVector_s {
        unsigned short distX, distY, distZ;

        DistVector_s( const unsigned short infty ) :
            distX( infty ), distY( infty ), distZ( infty ) {};

        DistVector_s( const unsigned short distX, const unsigned short distY,
            const unsigned short distZ ) :
            distX( distX ), distY( distY ), distZ( distZ ) {};

        DistVector_s operator + ( const DistVector_s & op2 ) {
            return DistVector_s( distX + op2.distX,
                distY + op2.distY, distZ + op2.distZ );
        }

        bool operator == ( const DistVector_s & op2 ) {
            return
                ( distX == op2.distX )
                && ( distY == op2.distY )
                && ( distZ == op2.distZ );
        }

        bool operator != ( const DistVector_s & op2 ) {
            return
            ( distX != op2.distX )
            || ( distY != op2.distY )
            || ( distZ != op2.distZ );
        }
    };

    static DistVector_s min( const DistVector_s & op1, const DistVector_s & op2 ) {

     if ( math::sqr( op1.distX ) + math::sqr( op1.distY ) + math::sqr( op1.distZ ) <=
          math::sqr( op2.distX ) + math::sqr( op2.distY ) + math::sqr( op2.distZ ) )
          return op1;
     else
          return op2;
    }

    class DistVectorField_t: public Volume_t<DistVector_s> {

    public:

        DistVectorField_t( int sizeX, int sizeY, int sizeZ,
            const DistVector_s & infty ) : Volume_t<DistVector_s>( sizeX,
            sizeY, sizeZ, infty ) {}

        /*void set( const int x, const int y, const int z,
              const DistVector_s & value ) {
            DistVector_s oval = this->get( x, y, z );
            if ( oval != value )
                std::cout << "( " << x << ", " << y << ", " << z << ") <- ("
                << value.distX << ", " << value.distY << ", " << value.distZ <<
                " )\n";
            Volume_t<DistVector_s>::set( x, y, z, value );
        }*/
    };

    typedef enum { ASC, DESC } Scandir_t;

    static void scanXLine( DistVectorField_t & dvField,
        int j, int k, const Scandir_t scandir );
    static void scanXYPlane( DistVectorField_t & dvField,
        int k, const Scandir_t scandir );
    static void scanVolume( DistVectorField_t & dvField,
        const Scandir_t scandir );
};

/** Class BitfieldReconstruction_t performs a volumetric reconstruction of
  * a solid using the modified Nooruddin/Turk (1999) method. The input
  * bitfield is taken as a point sampling of the boundary of the solid
  * with a defined density. */

class BitfieldReconstruction_t : public ScalarField_t<float> {

public :

    /**
     * Reconstruct a solid from a bitfield sampling its boundary. Delta
     * corresponds to the inverse of linear density. This means it should
     * be some measure of Euclidian distance between two points in the
     * sample, measured along the boundary.
     */
    BitfieldReconstruction_t( const Bitfield_t & from,
        double delta, double filterCutoffPeriod = 3.0 );

protected :

    struct Poll_s {
        unsigned char positives, negatives;

        bool operator == ( const Poll_s & op ) const {
            return ( positives == op.positives && negatives == op.negatives );
        }

        bool operator != ( const Poll_s & op ) const {
            return ( positives != op.positives || negatives != op.negatives );
        }

        Poll_s() : positives( 0 ), negatives( 0 ) {};
    };

    class VotingField_t: public Volume_t<Poll_s> {

    public:
        VotingField_t( const Bitfield_t & from )
            : Volume_t<Poll_s>( from.sizeX(), from.sizeY(), from.sizeZ(),
                Poll_s() ) {}
    };

    /**
     * process a single scanline, specified by a pair of iterators,
     * updating voting field along the way. This class does it via
     * a modified parity-count algorithm, based on intersections with delta
     * neighbourhood of boundary samples.
     */
    void scanline(
        const DistanceMap_t<double>::Giterator_t & begin,
        const DistanceMap_t<double>::Giterator_t & end,
        VotingField_t & vfield,
        const double delta );

    /**
     * determine the outcome of a poll. In this class, simple majority
     * wins.
     */
    float pollResult( const Poll_s & poll );
};

/* implementation follows */

/* class Volume_t<Value_t> */

template <typename Value_t>
Volume_t<Value_t>::Volume_t(
    const int sizeX, const int sizeY, const int sizeZ,
    const Value_t & initValue )
    : _sizeX( sizeX ), _sizeY( sizeY ), _sizeZ( sizeZ ),
      _initValue( initValue ) {

        _rootSize =  int( round( exp( log( 2.0 ) * ceil( log( std::max(
            sizeX, std::max( sizeY, sizeZ ) ) ) / log( 2.0 ) ) ) ) );

        _root = new Node_s( initValue );
}

template <typename Value_t>
Volume_t<Value_t>::~Volume_t() {
    delete _root;
}

template <typename Value_t>
Volume_t<Value_t>::Node_s::~Node_s() {

    if ( type == GRAY )
        for ( int i = 0; i < 8; i++ )
            delete subnodes[i];
}

template <typename Value_t>
Value_t Volume_t<Value_t>::get( int i, int j, int k ) const {

    if ( i < 0 || i >= _sizeX || j < 0 || j >= _sizeY || k < 0 || k >= _sizeZ )
        return _initValue;

    return _root->get( _rootSize, Position_s( i, j, k ) );
}

template <typename Value_t>
void Volume_t<Value_t>::set( int i, int j, int k, const Value_t & value ) {

    if ( i < 0 || i >= _sizeX || j < 0 || j >= _sizeY || k < 0 || k >= _sizeZ )
        return;

    _root->set( _rootSize, Position_s( i, j, k ), value );
}

/* class Volume_t<Value_t>::Giterator_t */

template <typename Value_t>
bool Volume_t<Value_t>::Giterator_t::operator < ( const Giterator_t & s ) {

    assert( diff == s.diff );
    Displacement_s df = s.pos - pos;
    if ( df.x * diff.x < 0 || df.y * diff.y < 0 || df.z * diff.z < 0 )
        return false;
    if ( abs( df.x > 0 ) || abs( df.y > 0 ) || abs( df.z ) > 0 )
        return true;
    else
        return false;
}

template <typename Value_t>
bool Volume_t<Value_t>::Giterator_t::operator <= ( const Giterator_t & s ) {

    assert( diff == s.diff );
    return ( *this < s ) || ( pos == s.pos );
}


template <typename Value_t>
typename Volume_t<Value_t>::Giterator_t Volume_t<Value_t>::gend(
    const Giterator_t & begin ) {

    float u = std::max( _sizeX, std::max( _sizeY, _sizeZ ) );
    float toss;

    // find closest clipping plane intersection
    // right
    if ( begin.diff.x > 0 ) {
        toss = ( _sizeX + 0.5 - begin.pos.x ) / begin.diff.x;
        if ( toss < u ) u = toss;
    }

    // left
    if ( begin.diff.x < 0 ) {
        toss = ( -1.5 - begin.pos.x ) / begin.diff.x;
        if ( toss < u ) u = toss;
    }

    // top
    if ( begin.diff.y > 0 ) {
         toss = ( _sizeY + 0.5 - begin.pos.y ) / begin.diff.y;
         if ( toss < u ) u = toss;
    }

    // bottom
    if ( begin.diff.y < 0 ) {
        toss = ( -1.5 - begin.pos.y ) / begin.diff.y;
        if ( toss < u ) u = toss;
    }

    // front
    if ( begin.diff.z > 0 ) {
        toss = ( _sizeZ + 0.5 - begin.pos.z ) / begin.diff.z;
        if ( toss < u ) u = toss;
    }

    // back
    if ( begin.diff.z < 0 ) {
        toss = ( -1.5 - begin.pos.z ) / begin.diff.z;
        if ( toss < u ) u = toss;
    }

    // done
    return Giterator_t( *this,
        Position_s(
            int( begin.pos.x + floor( u ) * begin.diff.x ),
            int( begin.pos.y + floor( u ) * begin.diff.y ),
            int( begin.pos.z + floor( u ) * begin.diff.z ) ),
        Displacement_s( begin.diff ) );
}

template <typename Value_t>
class std::set<typename Volume_t<Value_t>::Position_s>
    Volume_t<Value_t>::iteratorPositions( const Displacement_s & diff ) const {

    std::set<Volume_t<Value_t>::Position_s> retval;

    if ( diff.x != 0 ) {
        for ( int i = 0; i < _sizeY; i++ )
            for ( int j = 0; j < _sizeZ; j++ )
                retval.insert( Position_s(
                    diff.x > 0 ? 0 : _sizeX - 1, i, j ) );
    }

    if ( diff.y != 0 ) {
        for ( int i = 0; i < _sizeX; i++ )
            for ( int j = 0; j < _sizeZ; j++ )
                retval.insert( Position_s(
                    i, diff.y > 0 ? 0 : _sizeY - 1, j ) );
    }

    if ( diff.z != 0 ) {
        for ( int i = 0; i < _sizeX; i++ )
            for ( int j = 0; j < _sizeY; j++ )
                retval.insert( Position_s(
                    i, j, diff.z > 0 ? 0 : _sizeZ - 1 ) );
    }

    return retval;
}

/* class Volume_t<Value_t>::Node_s */

template <typename Value_t>
typename Volume_t<Value_t>::Node_s::Octant_t
    Volume_t<Value_t>::Node_s::findOctant(
        int nodeSize, const Position_s & pos ) {

    Octant_t retval = LBB;

    assert( pos.x < nodeSize && pos.y < nodeSize && pos.z < nodeSize );

    if ( pos.x >= ( nodeSize >> 1 ) ) retval = Octant_t( retval | OCT_X );
    if ( pos.y >= ( nodeSize >> 1 ) ) retval = Octant_t( retval | OCT_Y );
    if ( pos.z >= ( nodeSize >> 1 ) ) retval = Octant_t( retval | OCT_Z );

    return retval;
}

template <typename Value_t>
typename Volume_t<Value_t>::Position_s Volume_t<Value_t>::Node_s::toOctant(
    Octant_t octant, int nodeSize, const Position_s & pos ) {

    Position_s retval( pos );

    if ( octant & OCT_X ) retval.x -= ( nodeSize >> 1 );
    if ( octant & OCT_Y ) retval.y -= ( nodeSize >> 1 );
    if ( octant & OCT_Z ) retval.z -= ( nodeSize >> 1 );

    return retval;
}

template <typename Value_t>
typename Volume_t<Value_t>::Position_s Volume_t<Value_t>::Node_s::fromOctant(
    Octant_t octant, int nodeSize, const Position_s & pos ) {

    Position_s retval( pos );

    if ( octant & OCT_X ) retval.x += ( nodeSize >> 1 );
    if ( octant & OCT_Y ) retval.y += ( nodeSize >> 1 );
    if ( octant & OCT_Z ) retval.z += ( nodeSize >> 1 );

    return retval;
}

template <typename Value_t>
Value_t Volume_t<Value_t>::Node_s::get( int nodeSize, const Position_s & pos ) {

    if ( type == SOLID ) return value;

    Octant_t octant = findOctant( nodeSize, pos );

    return subnodes[ octant ]->get( nodeSize >> 1,
        toOctant( octant, nodeSize, pos ) );
}

template <typename Value_t>
void Volume_t<Value_t>::Node_s::set( int nodeSize, const Position_s & pos,
    const Value_t & value ) {

    if ( type == SOLID && this->value == value ) {
        return;
    }

    if ( type == SOLID && this->value != value ) {

        if ( nodeSize == 1 ) {

            this->value = value;
            return;

        } else {

            type = GRAY;
            for ( int i = 0; i < 8; i++ )
                subnodes[i] = new Node_s( this->value );
            // proceed to the code below

        }

    }

    if ( type == GRAY ) {

        // first set to the proper value
        Octant_t octant = findOctant( nodeSize, pos );
        subnodes[ octant ]->set( nodeSize >> 1,
                          toOctant( octant, nodeSize, pos ), value );

        // check for possible node collapse
        bool fullMatch = true;
        for ( int i = 0; i < 8; i++ )
            if ( subnodes[i]->type != SOLID || subnodes[i]->value != value )
                fullMatch = false;

        if ( fullMatch ) {

            type = SOLID; this->value = value;
            for ( int i = 0; i < 8; i++ ) delete subnodes[i];
        }

        return;
    }
}

/* class GeoVolume_t */

template <class Value_t>
GeoVolume_t<Value_t>::GeoVolume_t( const VolumeBase_t::FPosition_s & lower,
    const VolumeBase_t::FPosition_s & upper, const double voxelSize,
    const Value_t & initValue )
    : Volume_t<Value_t>(
        int( round( ( upper.x - lower.x ) / voxelSize ) ),
        int( round( ( upper.y - lower.y ) / voxelSize ) ),
        int( round( ( upper.z - lower.z ) / voxelSize ) ), initValue ),
     _lower( lower ), _upper( upper ), _voxelSize( voxelSize ) {

    // extents need to be modified to be voxelSize divisable
    _upper.x = _lower.x + this->_sizeX * _voxelSize;
    _upper.y = _lower.y + this->_sizeY * _voxelSize;
    _upper.z = _lower.z + this->_sizeZ * _voxelSize;
}

template <class Value_t>
Value_t GeoVolume_t<Value_t>::fget( const double x, const double y,
    const double z ) {

    typename VolumeBase_t::Position_s
        pos( geo2grid( VolumeBase_t::FPosition_s( x, y, z ) ) );

    return Volume_t<Value_t>::get( pos.x, pos.y, pos.z );
}

template <class Value_t>
void GeoVolume_t<Value_t>::fset( const double x, const double y,
    const double z, const Value_t & value ) {

    typename VolumeBase_t::Position_s
        pos( geo2grid( VolumeBase_t::FPosition_s( x, y, z ) ) );

    Volume_t<Value_t>::set( pos.x, pos.y, pos.z, value );
}

template <class Value_t>
typename VolumeBase_t::Position_s GeoVolume_t<Value_t>::geo2grid(
    const VolumeBase_t::FPosition_s & gpos ) {
    return typename Volume_t<Value_t>::Position_s(
        round( ( gpos.x - _lower.x ) / ( _upper.x - _lower.x )
            * this->_sizeX - 0.5 ),
        round( ( gpos.y - _lower.y ) / ( _upper.y - _lower.y )
            * this->_sizeY - 0.5 ),
        round( ( gpos.z - _lower.z ) / ( _upper.z - _lower.z )
            * this->_sizeZ - 0.5 ) );
}

template <class Value_t>
typename VolumeBase_t::FPosition_s GeoVolume_t<Value_t>::grid2geo(
    const typename VolumeBase_t::Position_s & pos ) {
    return VolumeBase_t::FPosition_s(
        _lower.x + ( pos.x + 0.5 ) / this->_sizeX * ( _upper.x - _lower.x ),
        _lower.y + ( pos.y + 0.5 ) / this->_sizeY * ( _upper.y - _lower.y ),
        _lower.z + ( pos.z + 0.5 ) / this->_sizeZ * ( _upper.z - _lower.z ) );
}

template <class Value_t>
typename VolumeBase_t::FPosition_s GeoVolume_t<Value_t>::grid2geo(
    const VolumeBase_t::FPosition_s & pos ) {
    return VolumeBase_t::FPosition_s(
        _lower.x + ( pos.x + 0.5 ) / this->_sizeX * ( _upper.x - _lower.x ),
        _lower.y + ( pos.y + 0.5 ) / this->_sizeY * ( _upper.y - _lower.y ),
        _lower.z + ( pos.z + 0.5 ) / this->_sizeZ * ( _upper.z - _lower.z ) );
}


/* class ScalarField_t */

template <class Value_t>
template <typename DstVolume_t>
    void ScalarField_t<Value_t>::filter( 
        const math::FIRFilter_t & filter,
        const VolumeBase_t::Displacement_s & diff,
        DstVolume_t & dstVolume ) {

        assert( this->sizeX() == dstVolume.sizeX() );
        assert( this->sizeY() == dstVolume.sizeY() );
        assert( this->sizeZ() == dstVolume.sizeZ() );
        assert( diff != VolumeBase_t::Displacement_s( 0, 0, 0 ) );

        std::set<VolumeBase_t::Position_s> poss = this->iteratorPositions( diff );

        BOOST_FOREACH( VolumeBase_t::Position_s pos, poss ) {

            typename ScalarField_t<Value_t>::Giterator_t sit( *this, pos, diff );
            typename ScalarField_t<Value_t>::Giterator_t send = gend( sit );
            typename DstVolume_t::Giterator_t dit( dstVolume, pos, diff );

            int rowSize = send - sit;

            for ( int x = 0; x < rowSize; x++ ) {
                dit.setValue( filter.convolute( sit, x, rowSize ) );
                ++sit; ++dit;
            }
        }
}

template <class Value_t>
    std::vector<typename GeoVolume_t<Value_t>::FPosition_s>
    ScalarField_t<Value_t>::getQuads( const Value_t & threshold,
        const SurfaceOrientation_t orientation ) {

    std::vector<typename GeoVolume_t<Value_t>::FPosition_s> retval;

    // iterate through all pixels
    for ( int i = 0; i < this->_sizeX; i++ )
        for ( int j = 0; j < this->_sizeY; j++ )
            for ( int k = 0; k < this->_sizeZ; k++ ) {

                // left
                if ( ( this->get( i, j, k ) > threshold
                     && this->get( i - 1, j, k ) <= threshold
                     && orientation == TO_MIN )
                     || ( this->get( i, j, k ) < threshold
                     && this->get ( i - 1, j, k ) >= threshold
                     && orientation == TO_MAX ) ) {

                        retval.push_back( grid2geo(
                            typename GeoVolume_t<Value_t>::FPosition_s(
                                i - 0.5, j - 0.5, k - 0.5 ) ) );
                        retval.push_back( grid2geo(
                            typename GeoVolume_t<Value_t>::FPosition_s(
                            i - 0.5, j - 0.5, k + 0.5 ) ) );
                        retval.push_back( grid2geo(
                            typename GeoVolume_t<Value_t>::FPosition_s(
                            i - 0.5, j + 0.5, k + 0.5 ) ) );
                        retval.push_back( grid2geo(
                            typename GeoVolume_t<Value_t>::FPosition_s(
                            i - 0.5, j + 0.5, k - 0.5 ) ) );
                }

                // right
                if ( ( this->get( i, j, k ) > threshold
                    && this->get( i + 1, j, k ) <= threshold
                    && orientation == TO_MIN )
                    || ( this->get( i, j, k ) < threshold
                    && this->get ( i + 1, j, k ) >= threshold
                    && orientation == TO_MAX ) ) {

                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j + 0.5, k - 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j + 0.5, k + 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j - 0.5, k + 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j - 0.5, k - 0.5 ) ) );
                    }

                //  bottom
                if ( ( this->get( i, j, k ) > threshold
                    && this->get( i, j - 1, k ) <= threshold
                    && orientation == TO_MIN )
                    || ( this->get( i, j, k ) < threshold
                    && this->get ( i, j - 1, k ) >= threshold
                    && orientation == TO_MAX ) ) {

                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i - 0.5, j - 0.5, k - 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j - 0.5, k - 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j - 0.5, k + 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i - 0.5, j - 0.5, k + 0.5 ) ) );
                }

                //  top
                if ( ( this->get( i, j, k ) > threshold
                    && this->get( i, j + 1, k ) <= threshold
                    && orientation == TO_MIN )
                    || ( this->get( i, j, k ) < threshold
                    && this->get ( i, j + 1, k ) >= threshold
                    && orientation == TO_MAX ) ) {

                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i - 0.5, j + 0.5, k + 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j + 0.5, k + 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j + 0.5, k - 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i - 0.5, j + 0.5, k - 0.5 ) ) );
                }

                // back
                if ( ( this->get( i, j, k ) > threshold
                    && this->get( i, j, k - 1 ) <= threshold
                    && orientation == TO_MIN )
                    || ( this->get( i, j, k ) < threshold
                    && this->get ( i, j, k - 1 ) >= threshold
                    && orientation == TO_MAX ) ) {

                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                        i - 0.5, j - 0.5, k - 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                        i - 0.5, j + 0.5, k - 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                        i + 0.5, j + 0.5, k - 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                        i + 0.5, j - 0.5, k - 0.5 ) ) );
                }

                // front
                if ( ( this->get( i, j, k ) > threshold
                    && this->get( i, j, k + 1 ) <= threshold
                    && orientation == TO_MIN )
                    || ( this->get( i, j, k ) < threshold
                    && this->get ( i, j, k + 1 ) >= threshold
                    && orientation == TO_MAX ) ) {

                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j - 0.5, k + 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i + 0.5, j + 0.5, k + 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i - 0.5, j + 0.5, k + 0.5 ) ) );
                    retval.push_back( grid2geo(
                        typename GeoVolume_t<Value_t>::FPosition_s(
                            i - 0.5, j - 0.5, k + 0.5 ) ) );
                }
            }

    return retval;
}

template <typename Value_t>
typename GeoVolume_t<Value_t>::FPosition_s ScalarField_t<Value_t>::interpolate(
    const typename GeoVolume_t<Value_t>::FPosition_s & p1,
    const Value_t & value1,
    const typename GeoVolume_t<Value_t>::FPosition_s & p2,
    const Value_t & value2,
    Value_t midval ) {


    double alpha = ( midval - value1 ) / ( value2 - value1 );

    return typename GeoVolume_t<Value_t>::FPosition_s(
        p1.x * ( 1.0 - alpha ) + p2.x * alpha,
        p1.y * ( 1.0 - alpha ) + p2.y * alpha,
        p1.z * ( 1.0 - alpha ) + p2.z * alpha );

}

template <typename Value_t>
void ScalarField_t<Value_t>::isoFromTetrahedron(
    std::vector<typename GeoVolume_t<Value_t>::FPosition_s> & retval,
    const typename GeoVolume_t<Value_t>::FPosition_s & vx0,
    const Value_t & value0,
    const typename GeoVolume_t<Value_t>::FPosition_s & vx1,
    const Value_t & value1,
    const typename GeoVolume_t<Value_t>::FPosition_s & vx2,
    const Value_t & value2,
    const typename GeoVolume_t<Value_t>::FPosition_s & vx3,
    const Value_t & value3,
    const Value_t & threshold, const SurfaceOrientation_t orientation ) {

    // case 0000, 1111
    if ( ( value0 > threshold && value1 > threshold && value2 > threshold
        && value3 > threshold )
        || ( value0 <= threshold && value1 <= threshold && value2 <= threshold
        && value3 <= threshold ) )
        return;

    // case 0001
    if ( ( ( value0 > threshold && value1 <= threshold && value2 <= threshold
        && value3 <= threshold && orientation == TO_MIN )
        || ( value0 <= threshold && value1 > threshold && value2 > threshold
        && value3 > threshold && orientation == TO_MAX ) ) ) {
        retval.push_back(
            interpolate( vx1, value1, vx0, value0, threshold ) );
        retval.push_back(
            interpolate( vx2, value2, vx0, value0, threshold ) );
        retval.push_back(
            interpolate( vx3, value3, vx0, value0, threshold ) );
        //std::cout << "0001\n";
    }

    // case 0010
    if ( ( ( value0 <= threshold && value1 > threshold && value2 <= threshold
        && value3 <= threshold && orientation == TO_MIN )
        || ( value0 > threshold && value1 <= threshold && value2 > threshold
        && value3 > threshold && orientation == TO_MAX ) ) ) {
        retval.push_back(
            interpolate( vx2, value2, vx1, value1, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx1, value1, threshold ) );
        retval.push_back(
            interpolate( vx3, value3, vx1, value1, threshold ) );
        //std::cout << "0010\n";
    }

    // case 0011
    if ( ( ( value0 > threshold && value1 > threshold && value2 <= threshold
        && value3 <= threshold && orientation == TO_MIN )
        || ( value0 <= threshold && value1 <= threshold && value2 > threshold
        && value3 > threshold && orientation == TO_MAX ) ) ) {
            retval.push_back(
                interpolate( vx1, value1, vx2, value2, threshold ) );
            retval.push_back(
                interpolate( vx0, value0, vx2, value2, threshold ) );
            retval.push_back(
                interpolate( vx1, value1, vx3, value3, threshold ) );
            retval.push_back(
                interpolate( vx1, value1, vx3, value3, threshold ) );
            retval.push_back(
                interpolate( vx0, value0, vx2, value2, threshold ) );
            retval.push_back(
                interpolate( vx0, value0, vx3, value3, threshold ) );
            //std::cout << "0011\n";
    }

    // case 0100
    if ( ( ( value0 <= threshold && value1 <= threshold && value2 > threshold
        && value3 <= threshold && orientation == TO_MIN )
        || ( value0 > threshold && value1 > threshold && value2 <= threshold
        && value3 > threshold && orientation == TO_MAX ) ) ) {
        retval.push_back(
            interpolate( vx1, value1, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx3, value3, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx2, value2, threshold ) );
        //std::cout << "0100\n";
    }

    // case 0101
    if ( ( ( value0 > threshold && value1 <= threshold && value2 > threshold
        && value3 <= threshold && orientation == TO_MIN )
        || ( value0 <= threshold && value1 > threshold && value2 <= threshold
        && value3 > threshold && orientation == TO_MAX ) ) ) {
        retval.push_back(
            interpolate( vx0, value0, vx1, value1, threshold ) );
        retval.push_back(
            interpolate( vx1, value1, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx2, value2, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx2, value2, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx1, value1, threshold ) );
        //std::cout << "0101\n";
    }

    // case 0110
    if ( ( ( value0 <= threshold && value1 > threshold && value2 > threshold
        && value3 <= threshold && orientation == TO_MIN )
        || ( value0 > threshold && value1 <= threshold && value2 <= threshold
        && value3 > threshold && orientation == TO_MAX ) ) ) {
        retval.push_back(
            interpolate( vx0, value0, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx1, value1, threshold ) );
        retval.push_back(
            interpolate( vx1, value1, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx1, value1, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx2, value2, vx3, value3, threshold ) );
        //std::cout << "0110\n";
    }

    // case 0111
    if ( ( ( value0 > threshold && value1 > threshold && value2 > threshold
        && value3 <= threshold && orientation == TO_MIN )
        || ( value0 <= threshold && value1 <= threshold && value2 <= threshold
        && value3 > threshold && orientation == TO_MAX ) ) ) {
        retval.push_back(
            interpolate( vx2, value2, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx1, value1, vx3, value3, threshold ) );
        //std::cout << "0111\n";
    }

    // case 1000
    if ( ( ( value0 <= threshold && value1 <= threshold && value2 <= threshold
        && value3 > threshold && orientation == TO_MIN )
        || ( value0 > threshold && value1 > threshold && value2 > threshold
        && value3 <= threshold && orientation == TO_MAX ) ) ) {
            retval.push_back(
                interpolate( vx2, value2, vx3, value3, threshold ) );
            retval.push_back(
                interpolate( vx1, value1, vx3, value3, threshold ) );
            retval.push_back(
                interpolate( vx0, value0, vx3, value3, threshold ) );
            //std::cout << "1000\n";
    }

    // case 1001
    if ( ( ( value0 > threshold && value1 <= threshold && value2 <= threshold
        && value3 > threshold && orientation == TO_MIN )
        || ( value0 <= threshold && value1 > threshold && value2 > threshold
        && value3 <= threshold && orientation == TO_MAX ) ) ) {
        retval.push_back(
            interpolate( vx0, value0, vx1, value1, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx1, value1, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx1, value1, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx2, value2, vx3, value3, threshold ) );
        //std::cout << "1001\n";
    }

    // case 1010
    if ( ( ( value0 <= threshold && value1 > threshold && value2 <= threshold
        && value3 > threshold && orientation == TO_MIN )
        || ( value0 >  threshold && value1 <= threshold && value2 > threshold
        && value3 <= threshold && orientation == TO_MAX ) ) ) {
        retval.push_back(
            interpolate( vx1, value1, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx1, value1, threshold ) );
        retval.push_back(
            interpolate( vx2, value2, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx2, value2, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx1, value1, threshold ) );
        //std::cout << "1010\n";
    }

    // case 1011
    if ( ( ( value0 > threshold && value1 > threshold && value2 <= threshold
        && value3 > threshold && orientation == TO_MIN )
        || ( value0 <= threshold && value1 <= threshold && value2 > threshold
        && value3 <= threshold && orientation == TO_MAX ) ) ) {

            retval.push_back(
                interpolate( vx3, value3, vx2, value2, threshold ) );
            retval.push_back(
                interpolate( vx1, value1, vx2, value2, threshold ) );
            retval.push_back(
                interpolate( vx0, value0, vx2, value2, threshold ) );
            //std::cout << "1011\n";
        }

    // case 1100
    if ( ( ( value0 <= threshold && value1 <= threshold && value2 > threshold
        && value3 > threshold && orientation == TO_MIN )
        || ( value0 > threshold && value1 > threshold && value2 <= threshold
        && value3 <= threshold && orientation == TO_MAX ) ) ) {
        retval.push_back(
            interpolate( vx0, value0, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx1, value1, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx1, value1, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx2, value2, threshold ) );
        retval.push_back(
            interpolate( vx1, value1, vx3, value3, threshold ) );
        retval.push_back(
            interpolate( vx0, value0, vx3, value3, threshold ) );
        //std::cout << "1100\n";
    }

    // case 1101
    if ( ( ( value0 > threshold && value1 <= threshold && value2 > threshold
        && value3 > threshold && orientation == TO_MIN )
        || ( value0 <= threshold && value1 > threshold && value2 <= threshold
        && value3 <= threshold && orientation == TO_MAX ) ) ) {
            retval.push_back(
                interpolate( vx0, value0, vx1, value1, threshold ) );
            retval.push_back(
                interpolate( vx2, value2, vx1, value1, threshold ) );
            retval.push_back(
                interpolate( vx3, value3, vx1, value1, threshold ) );
            //std::cout << "1101\n";
    }


    // case 1110
    if ( ( ( value0 <= threshold && value1 > threshold && value2 > threshold
        && value3 > threshold && orientation == TO_MIN )
        || ( value0 > threshold && value1 <= threshold && value2 <= threshold
        && value3 <= threshold && orientation == TO_MAX ) ) ) {
            retval.push_back(
                interpolate( vx2, value2, vx0, value0, threshold ) );
            retval.push_back(
                interpolate( vx1, value1, vx0, value0, threshold ) );
            retval.push_back(
                interpolate( vx3, value3, vx0, value0, threshold ) );
            //std::cout << "1110\n";
    }

}

template <typename Value_t>
std::vector<typename VolumeBase_t::FPosition_s>
ScalarField_t<Value_t>::isosurface( const Value_t & threshold,
            const SurfaceOrientation_t orientation ) {

    std::vector<typename GeoVolume_t<Value_t>::FPosition_s> retval;

    for ( int i = -1; i < this->_sizeX; i++ )
        for ( int j = -1; j < this->_sizeY; j++ )
            for ( int k = -1; k < this->_sizeZ; k++ ) {


                struct {
                    typename GeoVolume_t<Value_t>::FPosition_s vertex;
                    Value_t value;
                } vertexes[8];

                vertexes[0].vertex = grid2geo(
                    typename Volume_t<Value_t>::Position_s( i, j, k ) );
                vertexes[0].value = this->get( i, j, k );
                vertexes[1].vertex = grid2geo(
                    typename Volume_t<Value_t>::Position_s( i + 1, j, k ) );
                vertexes[1].value = this->get( i + 1, j, k );
                vertexes[2].vertex = grid2geo(
                    typename Volume_t<Value_t>::Position_s( i, j + 1, k ) );
                vertexes[2].value = this->get( i, j + 1, k );
                vertexes[3].vertex = grid2geo(
                    typename Volume_t<Value_t>::Position_s( i + 1, j + 1, k ) );
                vertexes[3].value = this->get( i + 1, j + 1, k );
                vertexes[4].vertex = grid2geo(
                    typename Volume_t<Value_t>::Position_s( i, j, k + 1 ) );
                vertexes[4].value = this->get( i, j, k + 1 );
                vertexes[5].vertex = grid2geo(
                    typename Volume_t<Value_t>::Position_s( i + 1, j, k + 1 ) );
                vertexes[5].value = this->get( i + 1, j, k + 1 );
                vertexes[6].vertex = grid2geo(
                    typename Volume_t<Value_t>::Position_s( i, j + 1, k + 1 ) );
                vertexes[6].value = this->get( i, j + 1, k + 1 );
                vertexes[7].vertex = grid2geo(
                    typename Volume_t<Value_t>::Position_s( i + 1, j + 1, k + 1 ) );
                vertexes[7].value = this->get( i + 1, j + 1, k + 1 );

                isoFromTetrahedron(
                    retval,
                    vertexes[0].vertex, vertexes[0].value,
                    vertexes[5].vertex, vertexes[5].value,
                    vertexes[7].vertex, vertexes[7].value,
                    vertexes[4].vertex, vertexes[4].value,
                    threshold, orientation );
                isoFromTetrahedron(
                    retval,
                    vertexes[0].vertex, vertexes[0].value,
                    vertexes[1].vertex, vertexes[1].value,
                    vertexes[7].vertex, vertexes[7].value,
                    vertexes[5].vertex, vertexes[5].value,
                    threshold, orientation );
                isoFromTetrahedron(
                    retval,
                    vertexes[0].vertex, vertexes[0].value,
                    vertexes[1].vertex, vertexes[1].value,
                    vertexes[3].vertex, vertexes[3].value,
                    vertexes[7].vertex, vertexes[7].value,
                    threshold, orientation );
                isoFromTetrahedron(
                    retval,
                    vertexes[0].vertex, vertexes[0].value,
                    vertexes[7].vertex, vertexes[7].value,
                    vertexes[6].vertex, vertexes[6].value,
                    vertexes[4].vertex, vertexes[4].value,
                    threshold, orientation );
                //if ( i == 0 && j == 0 && k == 1 )
                isoFromTetrahedron(
                    retval,
                    vertexes[0].vertex, vertexes[0].value,
                    vertexes[7].vertex, vertexes[7].value,
                    vertexes[2].vertex, vertexes[2].value,
                    vertexes[6].vertex, vertexes[6].value,
                    threshold, orientation );
                isoFromTetrahedron(
                    retval,
                    vertexes[0].vertex, vertexes[0].value,
                    vertexes[3].vertex, vertexes[3].value,
                    vertexes[2].vertex, vertexes[2].value,
                    vertexes[7].vertex, vertexes[7].value,
                    threshold, orientation );
            }

    return retval;
}


/** Class DistanceMap_t */

template <typename Value_t>
DistanceMap_t<Value_t>::DistanceMap_t( const Bitfield_t & bitfield,
    const Value_t initValue )
    : ScalarField_t<Value_t>( bitfield.lower(), bitfield.upper(),
      bitfield.voxelSize(), initValue ) {

    // initialize vector distance field (Danielsson's 4SED algorithm)
    DistVectorField_t dvField( this->_sizeX, this->_sizeY, this->_sizeZ,
        (unsigned short) ceil( initValue / this->_voxelSize ) );

    for ( int i = 0; i < dvField.sizeX(); i++ )
        for ( int j = 0; j < dvField.sizeY(); j++ )
            for ( int k = 0; k < dvField.sizeZ(); k++ )
                if ( bitfield.get( i, j, k ) )
                    dvField.set( i, j, k, DistVector_s( 0.0, 0.0, 0.0 ) );

    // perform scanning
    //std::cout << "Z ASC\n";
    for ( int k = 1; k < dvField.sizeZ(); k++ )
        scanXYPlane( dvField, k, ASC );
    //std::cout << "Z DESC\n";
    for ( int k = dvField.sizeZ() - 2; k >= 0; k-- )
        scanXYPlane( dvField, k, DESC );

    // compute distance map based on the vector field
    for ( int i = 0; i < dvField.sizeX(); i++ )
        for ( int j = 0; j < dvField.sizeY(); j++ )
            for ( int k = 0; k < dvField.sizeZ(); k++ ) {

                DistVector_s dv( dvField.get( i, j, k ) );
                Value_t dist = this->_voxelSize * sqrt(
                    math::sqr( dv.distX ) + math::sqr( dv.distY ) + math::sqr( dv.distZ ) );

                if ( dist < initValue )
                    set( i, j, k, dist );
            }

    // all done
}

template <typename Value_t>
void DistanceMap_t<Value_t>::scanXYPlane( DistVectorField_t & dvField,
    int k, const Scandir_t scandir ) {

    // z propagation
    for ( int i = 0; i < dvField.sizeX(); i++ )
        for ( int j = 0; j < dvField.sizeY(); j++ )
            if ( scandir == ASC ) {

                dvField.set( i, j, k, min(
                    dvField.get( i, j, k ),
                    dvField.get( i, j, k - 1 ) + DistVector_s( 0, 0, 1 ) ) );

            } else {

                dvField.set( i, j, k, min(
                    dvField.get( i, j, k ),
                    dvField.get( i, j, k + 1 ) + DistVector_s( 0, 0, 1 ) ) );
            }

    // xy propagation
    for ( int j = 1; j < dvField.sizeY(); j++ ) {
        //std::cout << "X ASC\n";
        scanXLine( dvField, j, k, ASC );
    }

    for ( int j = dvField.sizeY() - 2; j >= 0; j-- ) {
        //std::cout << "X DESC\n";
        scanXLine( dvField, j, k, DESC );
    }
}

template <typename Value_t>
void DistanceMap_t<Value_t>::scanXLine( DistVectorField_t & dvField,
    int j, int k, const Scandir_t scandir ) {

    // y propagation
    for ( int i = 0; i < dvField.sizeX(); i++ )
        if ( scandir == ASC )
            dvField.set( i, j, k, min(
                dvField.get( i, j, k ),
                dvField.get( i, j - 1, k ) + DistVector_s( 0, 1, 0 ) ) );
        else
            dvField.set( i, j, k, min(
                dvField.get( i, j, k ),
                dvField.get( i, j + 1, k ) + DistVector_s( 0, 1, 0 ) ) );

    // x propagation
    for ( int i = 1; i < dvField.sizeX(); i++ )
        dvField.set( i, j, k, min(
            dvField.get( i, j, k ),
            dvField.get( i - 1, j, k ) + DistVector_s( 1, 0, 0 ) ) );

    for ( int i = dvField.sizeX() - 2; i >= 0; i-- )
        dvField.set( i, j, k, min(
            dvField.get( i, j, k ),
            dvField.get( i + 1, j, k ) + DistVector_s( 1, 0, 0 ) ) );
}

} // namespace geometry

#endif

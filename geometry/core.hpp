/**
 * @file core.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Basic structures for geometric modeling - points and point containers
 */

#ifndef GEOMETRY_CORE_HPP
#define GEOMETRY_CORE_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <vector>

namespace ublas = boost::numeric::ublas;

namespace math {

class Point2 : public ublas::vector<double, ublas::bounded_array<double, 2> > {
public:
    Point2( const double x = 0.0, const double y = 0.0 )
        : ublas::vector<double, ublas::bounded_array<double, 2> >(2) {
        (*this)(0) = x; (*this)(1) = y;
    }

    template <class AE>
    Point2( const ublas::vector_expression<AE> & op )
        : ublas::vector<double, ublas::bounded_array<double, 2> >(2) {
        ublas::vector_assign<ublas::scalar_assign>(*this, op );
    }

    Point2( const ublas::vector<double> & op )
        : ublas::vector<double, ublas::bounded_array<double, 2> >(2) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
    }
};

class Point3 : public ublas::vector<double, ublas::bounded_array<double, 3> > {
public:
    Point3( const double x = 0.0, const double y = 0.0, const double z = 0.0 )
        : ublas::vector<double, ublas::bounded_array<double, 3> >(3) {
        (*this)(0) = x; (*this)(1) = y; (*this)(2) = z;
    }

    template <class AE>
    Point3( const ublas::vector_expression<AE> & op )
        : ublas::vector<double, ublas::bounded_array<double, 3> >(3) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
    }

    Point3( const ublas::vector<double> & op )
        : ublas::vector<double, ublas::bounded_array<double, 3> >(3) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
    }    
};

typedef std::vector<Point2> Points2;
typedef std::vector<Point3> Points3;


} // namespace math

#endif

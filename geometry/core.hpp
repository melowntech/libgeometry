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

namespace geometry {

class Point2 : public ublas::vector<double> {
public:
    Point2( const double x = 0.0, const double y = 0.0 ) : ublas::vector<double>(2) {
        (*this)(0) = x; (*this)(1) = y;
    }
};

class Point3 : public ublas::vector<double> {
public:
    Point3( const double x = 0.0, const double y = 0.0, const double z = 0.0 ) : ublas::vector<double>(3) {
        (*this)(0) = x; (*this)(1) = y; (*this)(2) = z;
    }
};

typedef std::vector<Point2> Points2;
typedef std::vector<Point3> Points3;


} // namespace geometry

#endif

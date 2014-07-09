/**
 *  @file smoothcloud.hpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Cloud smoothing stuff.
 *
 *  2014-07-09 (vasek) Stolen from window-mesh
 */

#ifndef geometry_smoothcloud_hpp_included_
#define geometry_smoothcloud_hpp_included_

#include <iostream>
#include <boost/program_options.hpp>

#include "math/geometry_core.hpp"

#include "./kdtree.hpp"

namespace geometry {

class CloudSmoother {
public:
    struct Params {
        /** Number of neihgbors to use for smoothing.
         */
        std::size_t neighbors;

        /** Radius of neighborhood to use for smoothing.
         */
        double radius;

        Params()
            : neighbors(10), radius(0.1)
        {}

        void configuration(const std::string &section
                           , boost::program_options::options_description
                           &config)
        {
            using boost::program_options::value;
            config.add_options()
                ((section + "neighbors").c_str(),
                 value(&neighbors)->default_value(neighbors),
                 "Number of point's neigbors to use.")
                ((section + "radius").c_str(),
                 value(&radius)->default_value(radius),
                 "Radius of point's neighborhood use.")
                ;
        }

        void configure(const boost::program_options::variables_map &) {}

        template <typename E, typename T>
        std::basic_ostream<E, T>& dump(std::basic_ostream<E, T> &os
                                       , const std::string &section
                                       = std::string()) const;

        bool feasible(std::size_t pointCount) const {
            return (neighbors <= pointCount) && (radius > 0.0);
        }
    };

    CloudSmoother(const Params &params = Params())
        : params_(params)
    {}

    const Params& params() const { return params_; }

    /** Smooth given cloud (in-place).
     */
    math::Points3 operator()(const math::Points3 &points) const;

    /** Smooth given cloud (in-place).
     *  Uses pre-generated kdtree
     */
    math::Points3 operator()(const math::Points3 &points
                             , const KdTree<math::Point3> &kdtree) const;

private:
    const Params params_;
};

// inline implementation

inline math::Points3 CloudSmoother::operator()(const math::Points3 &points)
    const
{
    if (!params_.feasible(points.size())) { return points; }
    return operator()(points, {points.begin(), points.end()});
}

template <typename E, typename T>
std::basic_ostream<E, T>&
CloudSmoother::Params::dump(std::basic_ostream<E, T> &os
                            , const std::string &section) const
{
    return os
        << section << "neighbors = " << neighbors << "\n"
        << section << "radius = " << radius << "\n"
        ;
}

} // namespace geometry

#endif //geometry_smoothcloud_hpp_included_

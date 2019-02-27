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

#include "kdtree.hpp"

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
            return neighbors && (neighbors <= pointCount) && (radius > 0.0);
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

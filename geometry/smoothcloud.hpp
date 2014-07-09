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
    return operator()(points, {points.begin(), points.end()});
}

} // namespace geometry

#endif //geometry_smoothcloud_hpp_included_

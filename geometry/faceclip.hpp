/**
 * @file faceclip.hpp
 * @author Jakub Cerveny <jakub.cerveny@ext.citationtech.net>
 *
 * Triangle clipping algorithm.
 */

#ifndef geometry_faceclip_hpp_included_
#define geometry_faceclip_hpp_included_

#include <vector>
#include <opencv2/core/core.hpp>
#include "math/geometry_core.hpp"

namespace geometry { namespace opencv {

//! Helper structure for clipping textured triangles.
struct ClipTriangle
{
    typedef cv::Point3d Point;
    typedef cv::Point2f TCoord;
    typedef std::vector<ClipTriangle> list;

    unsigned id1, id2;  // arbitrary user-specified IDs for the triangle
    Point pos[3]; // 3d position
    TCoord uv[3]; // texture coordinate

    ClipTriangle(unsigned id1, unsigned id2,
                 const Point &pos1, const Point &pos2, const Point &pos3,
                 const TCoord &uv1, const TCoord &uv2, const TCoord &uv3)
        : id1(id1), id2(id2)
    {
        pos[0] = pos1, pos[1] = pos2, pos[2] = pos3;
        uv[0] = uv1, uv[1] = uv2, uv[2] = uv3;
    }

    ClipTriangle(unsigned id1, unsigned id2,
                 const math::Point3 &pos1, const math::Point3 &pos2, const math::Point3 &pos3,
                 const math::Point2 &uv1, const math::Point2 &uv2, const math::Point2 &uv3)
        : id1(id1), id2(id2)
    {
        pos[0] = Point(pos1(0), pos1(1), pos1(2));
        pos[1] = Point(pos2(0), pos2(1), pos2(2));
        pos[2] = Point(pos3(0), pos3(1), pos3(2));

        uv[0] = TCoord(uv1(0), uv1(1));
        uv[1] = TCoord(uv2(0), uv2(1));
        uv[2] = TCoord(uv3(0), uv3(1));
    }
};

//! Clipping plane
struct ClipPlane
{
    cv::Point3d normal;
    double d;

    ClipPlane(double a, double b, double c, double d)
        : normal(a, b, c), d(d)
    {}

    ClipPlane()
        : normal(0, 0, 0), d(0)
    {}
};


//! Clips all triangles by a plane (i.e., removes parts on the negative side of
//! the plane), possibly producing some new triangles in the process.
//!
ClipTriangle::list clipTriangles(const ClipTriangle::list &triangles,
                                 const ClipPlane &plane);

} } // namespace geometry::opencv

#endif // geometry_faceclip_hpp_included_

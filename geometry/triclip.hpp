/**
 *  @file triclip.hpp
 *  @author Tomas Drinovsky <tomas.drinovsky@citationtech.net>
 *
 *  Triangle clipping stuff
 *  based on Jakub Cerveny's code
 */

#ifndef geometry_triclip_hpp_included_
#define geometry_triclip_hpp_included_


#include "dbglog/dbglog.hpp"

namespace geometry { 
/** Textured 3D triangle representation suitable for clipping algorithm.
 *
 * Triangle holds positions and tex. coordinates of its vertices.
 */
struct ClipTriangle
{
        typedef std::vector<ClipTriangle> list;
        // Indices of vertices.
        math::Point3 pos[3];
        // Indices of texture coordinates.
        math::Point2 uv[3];

        bool texCoordsAvailable;

        ClipTriangle(math::Point3 a, math::Point3 b
         , math::Point3 c)
        : texCoordsAvailable(false)
        {
            pos[0] = a, pos[1] = b, pos[2] = c;
        }

        ClipTriangle(math::Point3 a, math::Point3 b
         , math::Point3 c, math::Point2 ta
         , math::Point2 tb, math::Point2 tc)
        : texCoordsAvailable(true)
        {
            pos[0] = a, pos[1] = b, pos[2] = c;
            uv[0] = ta, uv[1] = tb, uv[2] = tc;
        }
};
/** Clipping plane
 *
 * Plane is defined by it's normal and the shift in the direction of the normal
 */
struct ClipPlane
{
        math::Point3 normal;
        double d;

        ClipPlane(double a, double b, double c, double d)
            : normal(a, b, c), d(d)
        {}

        ClipPlane()
            : normal(0, 0, 0), d(0)
        {}
};

namespace detail{
    inline double signedDistance(const math::Point3 &point, const ClipPlane &plane)
    {
        return math::dotProduct(point,plane.normal) - plane.d;
    }

    inline math::Point3 intersection(const math::Point3 &p1, const math::Point3 &p2,
                       const ClipPlane &plane, double& t)
    {
        // sort points to prevent numerical inaccuracies 
        math::Point3 sp1(std::max(p1,p2));
        math::Point3 sp2(std::min(p1,p2));

        double dot1 = math::dotProduct(sp1,plane.normal);
        double dot2 = math::dotProduct(sp2,plane.normal);
        double den = dot1 - dot2;

        // line parallel with plane, return the midpoint
        if (std::abs(den) < 1e-10) {
            t = 0.5;
            return (sp1 + sp2) * t;
        }

        t = (dot1 - plane.d) / den;
        return (1.0 - t)*sp1 + t*sp2;
    }
}

/** Clip triangles with given plane
 *
 * \param triangles list of triangles to be clipped
 * \param plane plane used to clip triangles 
 * \param triangleInfos additional triangle informations,
                        each newly created triangle will have same triangleInfo
                        as it's origin
 * \return list of the clipped triangles
 */
template<typename TriangleInfo>
ClipTriangle::list clipTriangles( const ClipTriangle::list &triangles
                                , const ClipPlane &plane
                                , std::vector<TriangleInfo> &triangleInfos)
{
    ClipTriangle::list result;
    std::vector<TriangleInfo> resultInfo;

    if (triangleInfos.size() && triangleInfos.size()!=triangles.size()) {
        LOGTHROW(err3, std::runtime_error)
            << "Triangle count and triangle informations count mismatch.";
    }

    for ( std::size_t tid=0; tid<triangles.size(); ++tid)
    {
        auto &tri(triangles[tid]);

        bool positive[3] = {
            detail::signedDistance(tri.pos[0], plane) >= 0,
            detail::signedDistance(tri.pos[1], plane) >= 0,
            detail::signedDistance(tri.pos[2], plane) >= 0
        };

        int count = 0;
        for (int i = 0; i < 3; i++) {
            if (positive[i]) count++;
        }

        // triangle completely on negative side - do nothing
        if (count == 0) continue;

        // trinagle completely on positive side - copy to result
        if (count == 3) {
            result.push_back(tri);
            if(triangleInfos.size()){
                resultInfo.emplace_back(triangleInfos[tid]);
            }
            continue;
        }

        int a, b, c;
        double t;

        // case 1: one vertex on positive side, just adjust the other two
        if (count == 1)
        {
            if (positive[0]) a = 0, b = 1, c = 2;
            else if (positive[1]) a = 1, b = 2, c = 0;
            else a = 2, b = 0, c = 1;

            auto x1pos(detail::intersection(tri.pos[a], tri.pos[b], plane, t));
            auto x1uv((1.0 - t)*tri.uv[a] + t*tri.uv[b]);

            auto x2pos(detail::intersection(tri.pos[c], tri.pos[a], plane, t));
            auto x2uv((1.0 - t)*tri.uv[c] + t*tri.uv[a]);

            result.emplace_back(tri.pos[a], x1pos, x2pos,
                                tri.uv[a],  x1uv,  x2uv);
            if(triangleInfos.size()){
                resultInfo.emplace_back(triangleInfos[tid]);
            }
        }
        // case 2: two vertices on positive side, adjust triangle and add one more
        else
        {
            if (!positive[0]) a = 0, b = 1, c = 2;
            else if (!positive[1]) a = 1, b = 2, c = 0;
            else a = 2, b = 0, c = 1;

            auto x1pos(detail::intersection(tri.pos[a], tri.pos[b], plane, t));
            auto x1uv((1.0 - t)*tri.uv[a] + t*tri.uv[b]);

            auto x2pos(detail::intersection(tri.pos[c], tri.pos[a], plane, t));
            auto x2uv((1.0 - t)*tri.uv[c] + t*tri.uv[a]);

            result.emplace_back(x1pos, tri.pos[b], tri.pos[c],
                                x1uv,  tri.uv[b],  tri.uv[c]);
            if(triangleInfos.size()){
                resultInfo.emplace_back(triangleInfos[tid]);
            }

            result.emplace_back(x1pos, tri.pos[c], x2pos,
                                x1uv,  tri.uv[c],  x2uv);
            if(triangleInfos.size()){
                resultInfo.emplace_back(triangleInfos[tid]);
            }
        }
    }

    triangleInfos.swap(resultInfo);
    return result;
}

} // namespace geometry
#endif // geometry_triclip_hpp_included_
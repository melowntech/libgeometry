/**
 *  @file geometry/parse-obj.hpp
 *  @author Vaclav Blazek <vaclav.blazek@ext.citationtech.net>
 *
 *  Boost.Spirit-based OBJ file format parser.
 */

#ifndef geometry_objparser_hpp_included_
#define geometry_objparser_hpp_included_

#include <iostream>
#include <vector>

#include <math/geometry_core.hpp>

namespace geometry {

class ObjParserBase {
public:
    struct Vector3d {
        double x, y, z;
        Vector3d() : x(), y(), z() {}

        operator math::Point3d() const {
            return math::Point3d(x, y, z);
        }
    };

    struct Facet {
        int va, vb, vc;
        int na, nb, nc;
        int ta, tb, tc;

        Facet() : va(), vb(), vc(), na(), nb(), nc(), ta(), tb(), tc() {}

        typedef std::vector<Facet> list;
    };

    virtual ~ObjParserBase() {}

    virtual void addVertex(const Vector3d&) = 0;
    virtual void addNormal(const Vector3d&) = 0;
    virtual void addFacet(const Facet&) = 0;

    bool parse(std::istream &is);
};

struct Obj : public ObjParserBase {
    typedef std::vector<math::Point3d> point_list;
    point_list vertices;
    point_list normals;
    Facet::list facets;

    void addVertex(const Vector3d &v) {
        vertices.push_back(v);
    }

    void addNormal(const Vector3d &n) {
        normals.push_back(n);
    }

    void addFacet(const Facet &f) {
        facets.push_back(f);
    }
};

} // namespace geometry

#endif // geometry_objparser_hpp_included_

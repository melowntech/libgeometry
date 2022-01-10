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
 *  @file geometry/parse-obj.hpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Boost.Spirit-based OBJ file format parser.
 */

#ifndef geometry_objparser_hpp_included_
#define geometry_objparser_hpp_included_

#include <iostream>
#include <vector>

#include <math/geometry_core.hpp>

#include <dbglog/dbglog.hpp>

namespace geometry {

class ObjParserBase {
public:
    struct Vector3d {
        double x, y, z;
        Vector3d() : x(), y(), z() {}

        operator math::Point3d() const {
            return math::Point3d(x, y, z);
        }

        bool operator<(const Vector3d &p) const;
    };

    struct Facet {
        int v[3];
        int n[3];
        int t[3];

        Facet() : v{0}, n{0}, t{0} {}

        typedef std::vector<Facet> list;
    };

    virtual ~ObjParserBase() {}

    virtual void addVertex(const Vector3d&) = 0;
    virtual void addTexture(const Vector3d&) = 0;
    virtual void addNormal(const Vector3d&) = 0;
    virtual void addFacet(const Facet&) = 0;

    virtual void materialLibrary(const std::string&) = 0;

    virtual void useMaterial(const std::string&) = 0;

    virtual void addObject(const std::string&) {};
    virtual void addGroup(const std::string&) {};

    // parses OBJ file from a stream, returns true if everything went OK
    bool parse(std::istream &is);

    // parses OBJ file from a stream, returns true if everything went OK
    bool parse(const boost::filesystem::path &path);
};

struct Obj : public ObjParserBase {
    math::Points3d vertices;
    math::Points3d texcoords;
    math::Points3d normals;
    Facet::list facets;

    void addVertex(const Vector3d &v) {
        vertices.push_back(v);
    }

    void addTexture(const Vector3d &t) {
        texcoords.push_back(t);
    }

    void addNormal(const Vector3d &n) {
        normals.push_back(n);
    }

    void addFacet(const Facet &f) {
        facets.push_back(f);
    }

    void materialLibrary(const std::string &l) {
        LOG(warn1) << "OBJ file contains material reference (mtllib " << l
                   << "); ignored by this simple reader.";
    }

    void useMaterial(const std::string &m) {
        LOG(warn1) << "OBJ file contains material reference (usemtl " << m
                   << "); ignored by this simple reader.";
    }
};

// io

template <typename E, typename T>
std::basic_ostream<E, T>&
operator<<(std::basic_ostream<E,T> &os
           , const geometry::ObjParserBase::Vector3d &v)
{
    return os << '[' << v.x << ' ' << v.y << ' ' << v.z << ']';
}

template <typename E, typename T>
std::basic_ostream<E, T>&
operator<<(std::basic_ostream<E,T> &os
           , const geometry::ObjParserBase::Facet &f)
{
    os << "f ";
    for (auto i : {0, 1, 2}) {
        if (i) os << ' ';
        if (f.v[i]) os << f.v[i];
        os << '/';
        if (f.t[i]) os << f.t[i];
        os << '/';
        if (f.n[i]) os << f.n[i];
    }

    return os;
}

// inlines
inline bool ObjParserBase::Vector3d::operator<(const Vector3d &p) const
{
    if (x < p.x) { return true; }
    if (p.x < x) { return false; }
    if (y < p.y) { return true; }
    if (p.y < y) { return false; }
    return (z < p.z);
}

} // namespace geometry

#endif // geometry_objparser_hpp_included_

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
 * @file math.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * 3D mesh representation and operations.
 */

#ifndef geometry_mesh_hpp_included_
#define geometry_mesh_hpp_included_

#include <memory>
#include <vector>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

#include "math/geometry.hpp"
#include "geometry/parse-obj.hpp"

namespace geometry {

/** Face in a mesh.
 */
struct Face {
    typedef std::vector<Face> list;

    unsigned int imageId;

    // Indices of vertices.
    math::Points3::size_type a, b, c;
    // Indices of texture coordinates.
    math::Points2::size_type ta, tb, tc;

    Face() : a(), b(), c(), ta(), tb(), tc() {}

    Face(math::Points3::size_type a, math::Points3::size_type b
         , math::Points3::size_type c)
        : imageId(), a(a), b(b), c(c), ta(), tb(), tc()
    {}

    Face(math::Points3::size_type a, math::Points3::size_type b
         , math::Points3::size_type c, math::Points2::size_type ta
         , math::Points2::size_type tb, math::Points2::size_type tc)
        : imageId(), a(a), b(b), c(c), ta(ta), tb(tb), tc(tc)
    {}

    /** Calculate normal of this face.
     */
    math::Point3 normal(const math::Points3 &vertices) const {
        return math::normalize
            (math::crossProduct(vertices[b] - vertices[a]
                                , vertices[c] - vertices[b]));
    }

    void clear() {
        a = b = c = 0;
        ta = tb = tc = 0;
    }

    /** Is this face not a triangle?
     */
    bool degenerate() const {
        return ((a == b) || (b == c) || (c == a));
    }
};

/** Textured 3D mesh representation.
 *
 * Vertices and texture coordinates are held in two arrays. Vertex at index I
 * has vertex coordinates at vertices[I] and texture coordinates at tCoords[I].
 *
 * Faces are defined as 3 indices to vertices and 3 indices to tCoords.
 */
struct Mesh {
    typedef std::shared_ptr<Mesh> pointer;

    /** vertices */
    math::Points3 vertices;

    /** per-vertex texture coordinates */
    math::Points2 tCoords;

    /** Faces (triplets of indices to vertices and texture coordinates) */
    Face::list faces;

    /** Face normal. */
    math::Point3 normal(const Face &face) const {
        return face.normal(vertices);
    }

    /** Add new face.
     */
    void addFace(math::Points3::size_type a, math::Points3::size_type b
                 , math::Points3::size_type c);

    void addFace(math::Points3::size_type a, math::Points3::size_type b
                 , math::Points3::size_type c, math::Points2::size_type ta
                 , math::Points2::size_type tb, math::Points2::size_type tc );

    /** First face point.
    */
    const math::Point3& a(const Face &face) const {
        return vertices[face.a];
    }

    /** Second face point.
     */
    const math::Point3& b(const Face &face) const {
        return vertices[face.b];
    }

    /** Third face point.
     */
    const math::Point3& c(const Face &face) const {
        return vertices[face.c];
    }

    /** Is given face not a triangle?
     */
    bool degenerate(const Face &face) const {
        return (face.degenerate() || (a(face) == b(face))
                || (b(face) == c(face)) || (c(face) == a(face)));
    }

    /** Are face indices ok?
     */
    bool good(const Face &face) const {
        return face.a < vertices.size() &&
               face.b < vertices.size() &&
               face.c < vertices.size();
    }

    /**
     * @brief provide mesh with skirt
     * @details skirt is a set quads pointing in the direction of the given
     * vector and attached to odd edges (edges adjacent to a single face).
     */
    void skirt( const math::Point3 & down = math::Point3( 0.0, 0.0, -1.0 ) );

    void sortFacesByImageId();

    /** Iterator that iteratex over face points (a->b->c->end)
     */
    struct FaceVertexConstIterator;

    FaceVertexConstIterator begin(const Face &face) const;

    FaceVertexConstIterator end(const Face&) const;

    /** Calculate face area.
     */
    double area(const Face &face);

    /** Calculate face barycenter.
     */
    math::Point3 barycenter(const Face &face);
};

inline void Mesh::addFace(math::Points3::size_type a
                          , math::Points3::size_type b
                          , math::Points3::size_type c)
{
    faces.emplace_back(a, b, c);
}

inline void Mesh::addFace(math::Points3::size_type a, math::Points3::size_type b
                 , math::Points3::size_type c, math::Points2::size_type ta
                 , math::Points2::size_type tb, math::Points2::size_type tc)
{
    faces.emplace_back(a, b, c, ta, tb, tc);
}


inline void Mesh::sortFacesByImageId()
{
    std::sort(faces.begin(), faces.end()
              , [](const Face &l, const Face &r) {
                  return l.imageId < r.imageId;
              });
}

struct Mesh::FaceVertexConstIterator {
    FaceVertexConstIterator() : face_(), vertices_(), index_() {}

    FaceVertexConstIterator(const Mesh &mesh, const Face &face)
        : face_(&face), vertices_(&mesh.vertices), index_()
    {}

    const math::Point3& operator*() const {
        switch (index_) {
        case 0: return (*vertices_)[face_->a];
        case 1: return (*vertices_)[face_->b];
        default: return (*vertices_)[face_->c];
        }
    }

    const math::Point3* operator->() const {
        switch (index_) {
        case 0: return &(*vertices_)[face_->a];
        case 1: return &(*vertices_)[face_->b];
        default: return &(*vertices_)[face_->c];
        }
    }

    FaceVertexConstIterator operator++() {
        ++index_;
        return *this;
    }

    bool operator==(const FaceVertexConstIterator &o) {
        if (((index_ > 2) || !face_) && !o.face_) {
            return true;
        }

        if (!face_ && ((o.index_ > 2) || !o.face_)) {
            return true;
        }

        return ((index_ == o.index_) && (face_ == o.face_)
                && (vertices_ == o.vertices_));
    }

    bool operator!=(const FaceVertexConstIterator &o) {
        return !operator==(o);
    }

private:
    const Face *face_;
    const math::Points3 *vertices_;
    unsigned int index_;
};

inline Mesh::FaceVertexConstIterator Mesh::begin(const Face &face) const {
    return FaceVertexConstIterator(*this, face);
}

inline Mesh::FaceVertexConstIterator Mesh::end(const Face&) const {
    return FaceVertexConstIterator();
}

} // namespace geometry
#endif // geometry_mesh_hpp_included_

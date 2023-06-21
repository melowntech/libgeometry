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
    typedef unsigned int index_type;

    unsigned int imageId;

    // Indices of vertices.
    index_type a, b, c;
    // Indices of texture coordinates.
    index_type ta, tb, tc;

    Face() : imageId(), a(), b(), c(), ta(), tb(), tc() {}

    Face(index_type a, index_type b, index_type c, unsigned int imageId = 0)
        : imageId(imageId), a(a), b(b), c(c), ta(), tb(), tc()
    {}

    Face(index_type a, index_type b, index_type c
         , index_type ta, index_type tb, index_type tc
         , unsigned int imageId = 0)
        : imageId(imageId), a(a), b(b), c(c), ta(ta), tb(tb), tc(tc)
    {}

    index_type vertex(const int idx) const {
        const index_type abc[] = { a, b, c };
        return abc[idx];
    }

    index_type& vertex(const int idx) {
        std::reference_wrapper<index_type> abc[] = { a, b, c };
        return abc[idx];
    }

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

    // needed by python bindings
    bool operator==(const Face &f) const;
    bool operator<(const Face &f) const;
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
    typedef std::vector<Mesh> list;

    /** vertices */
    std::vector <int> vertecesClass;
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
                 , math::Points3::size_type c, unsigned int imageId);

    void addFace(math::Points3::size_type a, math::Points3::size_type b
                 , math::Points3::size_type c, math::Points2::size_type ta
                 , math::Points2::size_type tb, math::Points2::size_type tc );

    void addFace(math::Points3::size_type a, math::Points3::size_type b
                 , math::Points3::size_type c, math::Points2::size_type ta
                 , math::Points2::size_type tb, math::Points2::size_type tc
                 , unsigned int imageId);

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

    /** First face texture point.
    */
    const math::Point2& ta(const Face &face) const {
        return tCoords[face.ta];
    }

    /** Second face texture point.
     */
    const math::Point2& tb(const Face &face) const {
        return tCoords[face.tb];
    }

    /** Third face texture point.
     */
    const math::Point2& tc(const Face &face) const {
        return tCoords[face.tc];
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

    /** Iterator that iterates over face points (a->b->c->end)
     */
    struct FaceVertexConstIterator;
    class FaceVertexConstIteratorRange;

    FaceVertexConstIterator begin(const Face &face) const;

    FaceVertexConstIterator end(const Face&) const;

    FaceVertexConstIteratorRange face(const Face &face) const;

    /** Calculate face area (in 3D space).
     */
    double area(const Face &face) const;

    /** Calculate face area (in UV space).;
     */
    double txArea(const Face &face) const;

    /** Calculate face barycenter.
     */
    math::Point3 barycenter(const Face &face) const;

    /** Calculate volume of a watertight mesh.
    */
    double volume() const;
};

// inlines

inline void Mesh::addFace(math::Points3::size_type a
                          , math::Points3::size_type b
                          , math::Points3::size_type c
                          , unsigned int imageId)
{
    faces.emplace_back(static_cast<geometry::Face::index_type>(a)
                      , static_cast<geometry::Face::index_type>(b)
                      , static_cast<geometry::Face::index_type>(c)
                      , imageId);
}

inline void Mesh::addFace(math::Points3::size_type a
                          , math::Points3::size_type b
                          , math::Points3::size_type c)
{
    faces.emplace_back(static_cast<geometry::Face::index_type>(a)
                      , static_cast<geometry::Face::index_type>(b)
                      , static_cast<geometry::Face::index_type>(c));
}

inline void
Mesh::addFace(math::Points3::size_type a, math::Points3::size_type b
              , math::Points3::size_type c, math::Points2::size_type ta
              , math::Points2::size_type tb, math::Points2::size_type tc)
{
    faces.emplace_back(static_cast<geometry::Face::index_type>(a)
                      , static_cast<geometry::Face::index_type>(b)
                      , static_cast<geometry::Face::index_type>(c)
                      , static_cast<geometry::Face::index_type>(ta)
                      , static_cast<geometry::Face::index_type>(tb)
                      , static_cast<geometry::Face::index_type>(tc));
}

inline void
Mesh::addFace(math::Points3::size_type a, math::Points3::size_type b
              , math::Points3::size_type c, math::Points2::size_type ta
              , math::Points2::size_type tb, math::Points2::size_type tc
              , unsigned int imageId)
{
    faces.emplace_back(static_cast<geometry::Face::index_type>(a)
                      , static_cast<geometry::Face::index_type>(b)
                      , static_cast<geometry::Face::index_type>(c)
                      , static_cast<geometry::Face::index_type>(ta)
                      , static_cast<geometry::Face::index_type>(tb)
                      , static_cast<geometry::Face::index_type>(tc)
                      , imageId);
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

class Mesh::FaceVertexConstIteratorRange {
public:
    FaceVertexConstIteratorRange(const Mesh &mesh, const Face &face)
        : begin_(mesh, face), end_()
    {}

    const FaceVertexConstIterator& begin() const { return begin_; }
    const FaceVertexConstIterator& end() const { return end_; }

private:
    FaceVertexConstIterator begin_;
    FaceVertexConstIterator end_;
};

inline Mesh::FaceVertexConstIteratorRange Mesh::face(const Face &face) const
{
    return { *this, face };
}

// inlines

inline bool Face::operator==(const Face &f) const {
    return
        (a == f.a) && (b == f.b) && (c == f.c)
        && (ta == f.ta) && (tb == f.tb) && (tc == f.tc)
        && (imageId == f.imageId)
        ;
}

inline bool Face::operator<(const Face &f) const
{
    if (a < f.a) { return true; }
    if (f.a < a) { return false; }
    if (b < f.b) { return true; }
    if (f.b < b) { return false; }
    if (c < f.c) { return true; }
    if (f.c < c) { return false; }

    if (ta < f.ta) { return true; }
    if (f.ta < ta) { return false; }
    if (tb < f.tb) { return true; }
    if (f.tb < tb) { return false; }
    if (tc < f.tc) { return true; }
    if (f.tc < tc) { return false; }

    return imageId < f.imageId;
}

} // namespace geometry

#endif // geometry_mesh_hpp_included_

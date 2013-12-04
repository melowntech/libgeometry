#ifndef mesh_hpp_included_
#define mesh_hpp_included_

#include <vector>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

#include "math/geometry.hpp"

namespace geometry {

struct Face {
    typedef std::vector<Face> list;

    unsigned int imageId;

    math::Points3::size_type a, b, c;
    math::Points2::size_type ta, tb, tc;

    Face() : a(), b(), c(), ta(), tb(), tc() {}

    Face(math::Points3::size_type a, math::Points3::size_type b
         , math::Points3::size_type c)
        : imageId(), a(a), b(b), c(c), ta(), tb(), tc()
    {}

    math::Point3 normal(const math::Points3 &vertices) const {
        return math::normalize
            (math::crossProduct(vertices[b] - vertices[a]
                                , vertices[c] - vertices[b]));
    }

    void clear() {
        a = b = c = 0;
        ta = tb = tc = 0;
    }

    bool degenerate() const {
        return ((a == b) || (b == c) || (c == a));
    }
};

struct Mesh {
    typedef boost::shared_ptr<Mesh> pointer;

    /** vertices */
    math::Points3 vertices;

    /** per-vertex texture coordinates */
    math::Points2 tCoords;

    /** Faces (triplets of indices to vertices and texture coordinates) */
    Face::list faces;

    math::Point3 normal(const Face &face) const {
        return face.normal(vertices);
    }

    void addFace(math::Points3::size_type a, math::Points3::size_type b
                 , math::Points3::size_type c);

    const math::Point3& a(const Face &face) const {
        return vertices[face.a];
    }

    const math::Point3& b(const Face &face) const {
        return vertices[face.b];
    }

    const math::Point3& c(const Face &face) const {
        return vertices[face.c];
    }

    bool degenerate(const Face &face) const {
        return (face.degenerate() || (a(face) == b(face))
                || (b(face) == c(face)) || (c(face) == a(face)));
    }

    bool good(const Face &face) const {
        return face.a < vertices.size() &&
               face.b < vertices.size() &&
               face.c < vertices.size();
    }

    void sortFacesByImageId();

    struct FaceVertexConstIterator;

    FaceVertexConstIterator begin(const Face &face) const;

    FaceVertexConstIterator end(const Face&) const;
};

inline void Mesh::addFace(math::Points3::size_type a
                          , math::Points3::size_type b
                          , math::Points3::size_type c)
{
    faces.emplace_back(a, b, c);
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
#endif // mesh_hpp_included_

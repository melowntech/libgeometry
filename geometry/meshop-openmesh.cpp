/**
 * @file math.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * 3D mesh operations
 */

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>

#include "./meshop.hpp"

namespace geometry {

/* code moved from the window-mesh utility */

namespace {

typedef OpenMesh::TriMesh_ArrayKernelT<> OMMesh;
typedef OpenMesh::Decimater::DecimaterT<OMMesh> Decimator;
typedef OpenMesh::Decimater::ModQuadricT<OMMesh>::Handle HModQuadric;

void toOpenMesh(const geometry::Mesh &mesh, OMMesh& omMesh)
{
    // create OpenMesh vertices
    std::vector<OMMesh::VertexHandle> handles;
    handles.reserve(mesh.vertices.size());
    for (const auto& v : mesh.vertices)
        handles.emplace_back(
            omMesh.add_vertex(OMMesh::Point(v(0), v(1), v(2))) );

    // create OpenMesh faces
    for (const geometry::Face& face : mesh.faces)
    {
        OMMesh::VertexHandle face_handles[3] = {
            handles[face.a], handles[face.b], handles[face.c]
        };
        omMesh.add_face(face_handles, 3);
    }
}

void fromOpenMesh(const OMMesh& omMesh, geometry::Mesh& mesh)
{
    // create our vertices
    for (auto v_it = omMesh.vertices_begin();
              v_it != omMesh.vertices_end();
              ++v_it)
    {
        const auto& pt = omMesh.point(v_it.handle());
        mesh.vertices.emplace_back(pt[0], pt[1], pt[2]);
    }

    // create our faces
    for (auto f_it = omMesh.faces_begin();
              f_it != omMesh.faces_end();
              ++f_it)
    {
        auto fv_it(omMesh.cfv_iter(f_it.handle()));
        int index[3], *p = index;
        for ( ; fv_it; ++fv_it) {
            *p++ = fv_it.handle().idx();
        }
        mesh.faces.emplace_back(index[0], index[1], index[2]);
    }
}

void lockCorners(OMMesh &omMesh)
{
    // get bounding box
    double box[2][2] = { {INFINITY, -INFINITY}, {INFINITY, -INFINITY} };

    for (auto v_it = omMesh.vertices_begin();
              v_it != omMesh.vertices_end();  ++v_it)
    {
        const auto& pt = omMesh.point(v_it.handle());
        for (int i = 0; i < 2; i++) {
            if (pt[i] < box[i][0]) box[i][0] = pt[i];
            if (pt[i] > box[i][1]) box[i][1] = pt[i];
        }
    }

    // find vertices nearest to the four corners of the bounding box
    struct { double dist; OMMesh::VertexHandle handle; } corners[2][2]
        = {{{INFINITY, {}}, {INFINITY, {}}}, {{INFINITY, {}}, {INFINITY, {}}}};

    omMesh.request_vertex_status();
    for (auto v_it = omMesh.vertices_begin();
              v_it != omMesh.vertices_end();  ++v_it)
    {
        const auto& pt = omMesh.point(v_it.handle());

        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            double dist = std::hypot(pt[0] - box[0][i], pt[1] - box[1][j]);
            if (dist < corners[i][j].dist) {
                corners[i][j].dist = dist;
                corners[i][j].handle = v_it.handle();
            }
        }
    }

    // lock the corners
    for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
    {
        auto handle = corners[i][j].handle;
        if (handle == OMMesh::VertexHandle()) continue;
        LOG(info2) << "Locking corner vertex " << handle.idx();
        omMesh.status(handle).set_locked(true);
    }
}

} // namespace

Mesh::pointer simplify(const Mesh &mesh, int faceCount)
{
    OMMesh omMesh;
    toOpenMesh(mesh, omMesh);

    // lock the corner vertices of the window to prevent simplifying the corners
    lockCorners(omMesh);

    Decimator decimator(omMesh);
    HModQuadric hModQuadric; // collapse priority based on vertex error quadric
    decimator.add(hModQuadric);
    decimator.initialize();

    decimator.decimate_to_faces(0, faceCount);
    omMesh.garbage_collection();

    auto newMesh(std::make_shared<geometry::Mesh>());
    fromOpenMesh(omMesh, *newMesh);
    return newMesh;
}

} // namespace geometry

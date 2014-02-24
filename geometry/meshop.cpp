/**
 * @file math.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * 3D mesh operations
 */

#include "./meshop.hpp"

namespace geometry {

Obj asObj(const Mesh &mesh)
{
    Obj obj;
    for ( math::Point3 vertex : mesh.vertices ) {

        geometry::Obj::Vector3d overtex;
        overtex.x = vertex[0]; overtex.y = vertex[1]; overtex.z = vertex[2];

        obj.addVertex( overtex );

        //LOG( info1 ) << "[" << overtex.x << ", " << overtex.y << ", "
        //    << overtex.z << "]";
    }

    for ( math::Point2 texture : mesh.tCoords ) {

        geometry::Obj::Vector3d otexture;
        otexture.x = texture[0]; otexture.y = texture[1]; otexture.z = 0.0;

        obj.addTexture( otexture );
    }

    for ( geometry::Face face : mesh.faces ) {

        geometry::Obj::Facet facet;

        facet.v[0] = face.a; facet.v[1] = face.b; facet.v[2] = face.c;
        facet.t[0] = face.ta; facet.t[1] = face.tb; facet.t[2] = face.tc;

        obj.addFacet( facet );
    }

    return obj;
}

void saveAsObj(const Mesh &mesh, const boost::filesystem::path &filepath
               , const std::string &mtlName)
{
    LOG(info4) << "Saving mesh to file <" << filepath << ">.";

    std::ofstream out(filepath.string().c_str());
    out.setf(std::ios::scientific, std::ios::floatfield);

    out << "mtllib " << mtlName << '\n';

    for (const auto &vertex : mesh.vertices) {
        out << "v " << vertex(0) << ' ' << vertex(1) << ' '  << vertex(2)
            << '\n';
    }

    for (const auto &tCoord : mesh.tCoords) {
        out << "vt " << tCoord(0) << ' ' << tCoord(1) << '\n';
    }

    unsigned int currentImageId(static_cast<unsigned int>(-1));

    for (const auto &face : mesh.faces) {
        if (face.degenerate()) {
            continue;
        }
        if (face.imageId != currentImageId) {
            out << "usemtl " << face.imageId << '\n';
            currentImageId = face.imageId;
        }

        out << "f " << face.a + 1 << '/' << face.ta + 1 << "/ "
            << face.b + 1 << '/' << face.tb + 1 << "/ "
            << face.c + 1 << '/' << face.tc + 1 << "/\n";
    }

    if (!out) {
        LOGTHROW(err3, std::runtime_error)
            << "Unable to save mesh to <" << filepath << ">.";
    }
}

} // namespace geometry

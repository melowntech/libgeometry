/**
 * @file math.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * 3D mesh operations
 */

#include "./meshop.hpp"
#include "./parse-obj.hpp"

#include "utility/expect.hpp"

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

Mesh loadPly( const boost::filesystem::path &filename )
{
    std::ifstream f(filename.native());
    if (!f.good()) {
        LOGTHROW(err2, std::runtime_error)
                << "Can't open " << filename;
    }

    f.exceptions(std::ios::badbit | std::ios::failbit);

    // read header
    std::string line;
    int nvert = -1, ntris = -1;
    do {
        if (getline(f, line).eof()) break;
        sscanf(line.c_str(), "element vertex %d", &nvert);
        sscanf(line.c_str(), "element face %d", &ntris);
    } while (line != "end_header");

    if (nvert < 0 || ntris < 0) {
        LOGTHROW(err2, std::runtime_error)
                << filename << ": unknown PLY format.";
    }

    Mesh mesh;

    // load points
    for (int i = 0; i < nvert; i++) {
        double x, y, z;
        f >> x >> y >> z;
        mesh.vertices.emplace_back(x, y, z);
    }

    // load triangles
    for (int i = 0; i < ntris; i++) {
        int n, a, b, c;
        f >> n >> a >> b >> c;
        utility::expect(n == 3, "Only triangles are supported in PLY files.");
        mesh.faces.emplace_back(a, b, c);
    }

    return mesh;
}

Mesh loadObj( const boost::filesystem::path &filename )
{
    struct Obj2MeshParser : public ObjParserBase {
        void addVertex( const Vector3d &v ) {
            mesh.vertices.push_back(v);
        }

        void addTexture( const Vector3d &t ) {
            mesh.tCoords.emplace_back(t.x, t.y);
        }

        void addFacet( const Facet &f ) {
            mesh.addFace( f.v[0], f.v[1], f.v[2]
                        , f.t[0], f.t[1], f.t[2] );
        }

        void addNormal( const Vector3d& ) { }
        void materialLibrary(const std::string&) { }
        void useMaterial(const std::string&) { }

        Mesh mesh;
    } parser_;

    std::ifstream file;
    file.exceptions(std::ios::badbit | std::ios::failbit);
    file.open(filename.string());

    parser_.parse(file);

    return parser_.mesh;
}

} // namespace geometry

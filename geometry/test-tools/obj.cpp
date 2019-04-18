#include <cstdlib>

#include "dbglog/dbglog.hpp"

#include "geometry/parse-obj.hpp"

int main(int argc, char *argv[])
{
    if (argc != 2) {
        std::cerr << "Invalid number of arguments." << std::endl;
        return EXIT_FAILURE;
    }

    struct Loader : geometry::ObjParserBase {
        virtual void addVertex(const Vector3d &v) {
            LOG(info4) << "vertex: (" << v.x << ", " << v.y << ", "
                       << v.z << ")";
        }
        virtual void addTexture(const Vector3d &t) {
            LOG(info4) << "texture: (" << t.x << ", " << t.y << ", "
                       << t.z << ")";
        }
        virtual void addNormal(const Vector3d &n) {
            LOG(info4) << "normal: (" << n.x << ", " << n.y << ", "
                       << n.z << ")";
        }
        virtual void addFacet(const Facet &f) {
            LOG(info4) << "facet: v=(" << f.v[0] << ", " << f.v[1] << ", "
                       << f.v[2] << "), t=("
                       << "facet: t=(" << f.t[0] << ", " << f.t[1] << ", "
                       << f.t[2] << "), b=("
                       << "facet: n=(" << f.n[0] << ", " << f.n[1] << ", "
                       << f.n[2] << ")";
        }

        virtual void materialLibrary(const std::string &path) {
            LOG(info4) << "material library: " << path;
        }

        virtual void useMaterial(const std::string &name) {
            LOG(info4) << "material: " << name;
        }
    };

    Loader loader;
    loader.parse(argv[1]);

    return EXIT_SUCCESS;
}

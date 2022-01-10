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
        void addVertex(const Vector3d &v) override {
            LOG(info4) << "vertex: (" << v.x << ", " << v.y << ", "
                       << v.z << ")";
        }
        void addTexture(const Vector3d &t) override {
            LOG(info4) << "texture: (" << t.x << ", " << t.y << ", "
                       << t.z << ")";
        }
        void addNormal(const Vector3d &n) override {
            LOG(info4) << "normal: (" << n.x << ", " << n.y << ", "
                       << n.z << ")";
        }
        void addFacet(const Facet &f) override {
            LOG(info4) << "facet: v=(" << f.v[0] << ", " << f.v[1] << ", "
                       << f.v[2] << "), t=("
                       << "facet: t=(" << f.t[0] << ", " << f.t[1] << ", "
                       << f.t[2] << "), b=("
                       << "facet: n=(" << f.n[0] << ", " << f.n[1] << ", "
                       << f.n[2] << ")";
        }

        void materialLibrary(const std::string &path) override {
            LOG(info4) << "material library: " << path;
        }

        void useMaterial(const std::string &name) override{
            LOG(info4) << "material: " << name;
        }

        void addObject(const std::string &name) override {
            LOG(info4) << "object: " << name;
        }

        void addGroup(const std::string &name) override {
            LOG(info4) << "group: " << name;
        }
    };

    Loader loader;
    auto res(loader.parse(argv[1]));
    if (!res) {
        LOG(fatal) << "Failed to parse input file.";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

#include <stdexcept>
#include <fstream>

#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <dbglog/dbglog.hpp>

#include "mesh.hpp"

namespace fs = boost::filesystem;

namespace math {

template<class Archive>
inline void serialize(Archive &ar, Point3 &p, const unsigned int version)
{
    (void) version;

    ar & p(0) & p(1) & p(2);
}

template<class Archive>
inline void serialize(Archive &ar, Point2 &p, const unsigned int version)
{
    (void) version;

    ar & p(0) & p(1);
}

} // namespace math

namespace geometry {

const unsigned int DATA_DUMP_VERSION(1);

void Mesh::saveObj(const boost::filesystem::path &filepath
                   , const std::string &mtlName) const
{
    LOG(info4) << "Saving mesh to file <" << filepath << ">.";

    std::ofstream out(filepath.string().c_str());
    out.setf(std::ios::scientific, std::ios::floatfield);

    out << "mtllib " << mtlName << '\n';

    for (const auto &vertex : vertices) {
        out << "v " << vertex(0) << ' ' << vertex(1) << ' '  << vertex(2)
            << '\n';
    }

    for (const auto &tCoord : tCoords) {
        out << "vt " << tCoord(0) << ' ' << tCoord(1) << '\n';
    }

    unsigned int currentImageId(static_cast<unsigned int>(-1));

    for (const auto &face : faces) {
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

template<class Archive>
inline void serialize(Archive &ar, Face &f, const unsigned int version)
{
    (void) version;

    ar & f.imageId & f.a & f.b & f.c & f.ta & f.tb & f.tc;
}


void save(const fs::path &path, const std::vector<std::string> &imagePaths
          , const Mesh &mesh)
{
    std::ofstream ofs(path.string().c_str(), std::ios::binary);
    if (ofs.fail()) {
        // TODO: use custom serialization exception
        LOGTHROW(err1, std::runtime_error)
            << "Cannot open mesh data file " << path << '.';
    }

    boost::archive::binary_oarchive oa(ofs);
    oa & DATA_DUMP_VERSION & imagePaths & mesh.vertices
        & mesh.tCoords & mesh.faces;
}

void load(const fs::path &path, std::vector<std::string> &imagePaths
          , Mesh &mesh)
{
    std::ifstream ifs(path.string().c_str(), std::ios::binary);
    if (ifs.fail()) {
        // TODO: use custom serialization exception
        LOGTHROW(err1, std::runtime_error)
            << "Cannot open mehs data file " << path << '.';
    }

    boost::archive::binary_iarchive ia(ifs);
    unsigned int version;
    ia & version;
    if (DATA_DUMP_VERSION != version) {
        // TODO: use custom serialization exception
        LOGTHROW(err1, std::runtime_error)
            << "Wrong mesh data version: "
            << " expected " << DATA_DUMP_VERSION << ", got " << version
            << " (file: " << path << ").";
    }

    imagePaths.clear();
    mesh.vertices.clear();
    mesh.tCoords.clear();
    mesh.faces.clear();

    ia & imagePaths & mesh.vertices & mesh.tCoords & mesh.faces;
}

} //namespace geometry

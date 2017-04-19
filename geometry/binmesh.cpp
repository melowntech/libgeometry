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
#include "utility/binaryio.hpp"

#include "./error.hpp"
#include "./binmesh.hpp"

namespace bin = utility::binaryio;

namespace geometry {

namespace {

const char MAGIC[8] = {'B', 'I', 'N', '.', 'M', 'E', 'S', 'H'};
const unsigned VERSION = 1;

} // namespace

void writeBinaryMesh(std::ostream &f, const geometry::Obj &mesh)
{
    math::Extents3 bbox(INFINITY, INFINITY, INFINITY,
                        -INFINITY, -INFINITY, -INFINITY);
    for (const auto &v : mesh.vertices) {
        update(bbox, v);
    }

    // write header
    bin::write(f, MAGIC);
    bin::write(f, uint32_t(VERSION));

    for (int i = 0; i < 3; i++) {
        bin::write(f, double(bbox.ll(i)));
    }
    for (int i = 0; i < 3; i++) {
        bin::write(f, double(bbox.ur(i)));
    }

    // write vertices
    int nv(std::min(mesh.vertices.size(), size_t(USHRT_MAX)));
    bin::write(f, uint16_t(nv));

    for (int i = 0; i < nv; i++) {
        const auto &v(mesh.vertices[i]);
        for (int j = 0; j < 3; j++) {
            double coord = (v(j) - bbox.ll(j)) / (bbox.ur(j) - bbox.ll(j));
            bin::write(f, uint16_t(round(coord * USHRT_MAX)));
        }
    }

    // write texture coords
    int ntv(std::min(mesh.texcoords.size(), size_t(USHRT_MAX)));
    bin::write(f, uint16_t(ntv));

    for (int i = 0; i < ntv; i++) {
        const auto &t(mesh.texcoords[i]);
        for (int j = 0; j < 2; j++) {
            bin::write(f, uint16_t(round(t(j) * USHRT_MAX)));
        }
    }

    // write faces
    int nf(std::min(mesh.facets.size(), size_t(USHRT_MAX)));
    bin::write(f, uint16_t(nf));

    for (int i = 0; i < nf; i++) {
        const auto &face(mesh.facets[i]);
        for (int j = 0; j < 3; j++) {
            bin::write(f, uint16_t(face.v[j]));
        }
        for (int j = 0; j < 3; j++) {
            bin::write(f, uint16_t(face.t[j]));
        }
    }
}

void writeBinaryMesh(const boost::filesystem::path &path,
                     const geometry::Obj &mesh)
{
    utility::ofstreambuf f(path.string());
    writeBinaryMesh(f, mesh);
    f.close();
}

geometry::Obj loadBinaryMesh(std::istream &f
                             , const boost::filesystem::path &path
                             , MeshInfo *meshInfo)
{
    // load header and check version
    char magic[8];
    uint32_t version;

    bin::read(f, magic);
    bin::read(f, version);

    if (std::memcmp(magic, MAGIC, sizeof(MAGIC))) {
        LOGTHROW(err1, BadFileFormat)
                << "File " << path << " is not a binary mesh file.";
    }
    if (version > VERSION) {
        LOGTHROW(err1, VersionError) << "File " << path
                << " has unsupported version (" << version << ").";
    }

    // load bounding box
    math::Extents3 bbox;
    for (int i = 0; i < 3; i++) {
        bin::read(f, bbox.ll(i));
    }
    for (int i = 0; i < 3; i++) {
        bin::read(f, bbox.ur(i));
    }

    geometry::Obj mesh;
    uint16_t count, coord, index;

    // load vertices
    bin::read(f, count);
    mesh.vertices.reserve(count);
    for (int k = 0; k < count; k++) {
        math::Point3 vertex;
        for (int i = 0; i < 3; i++) {
            bin::read(f, coord);
            double c(coord);
            c /= USHRT_MAX;
            vertex(i) = bbox.ll(i) + c * (bbox.ur(i) - bbox.ll(i));
        }
        mesh.vertices.push_back(vertex);
    }

    // load texcoords
    bin::read(f, count);
    mesh.texcoords.reserve(count);
    for (int k = 0; k < count; k++) {
        math::Point2 tex;
        for (int i = 0; i < 2; i++) {
            bin::read(f, coord);
            tex(i) = double(coord) / USHRT_MAX;
        }
        mesh.texcoords.emplace_back(tex(0), tex(1), 0.0);
    }

    // load faces
    bin::read(f, count);
    mesh.facets.reserve(count);
    for (int k = 0; k < count; k++) {
        geometry::Obj::Facet face;
        for (int i = 0; i < 3; i++) {
            bin::read(f, index);
            face.v[i] = index;
        }
        for (int i = 0; i < 3; i++) {
            bin::read(f, index);
            face.t[i] = index;
        }
        mesh.facets.push_back(face);
    }

    if (meshInfo) {
        meshInfo->bbox = bbox;
        meshInfo->vertexCount = mesh.vertices.size();
        meshInfo->faceCount = mesh.facets.size();
        meshInfo->texCoordCount = mesh.texcoords.size();
    }

    return mesh;
}

geometry::Obj loadBinaryMesh(const boost::filesystem::path &path
                             , MeshInfo *meshInfo)
{
    utility::ifstreambuf f(path.string());
    auto mesh(loadBinaryMesh(f, path, meshInfo));
    f.close();
    return mesh;
}

} // namespace geometry

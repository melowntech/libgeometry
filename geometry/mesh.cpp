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
#include <stdexcept>
#include <fstream>

#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <dbglog/dbglog.hpp>

#include "mesh.hpp"
#include "utility/openmp.hpp"
#include "utility/expect.hpp"

namespace fs = boost::filesystem;


namespace geometry {


void Mesh::computeBoundary() {
    // Face-face table is required for computing mesh boundary.
    utility::expect(connectivity.faceFaceTable.size() == faces.size());

    boundary.clear();
    boundary.resize(vertices.size(), false);

    for (uint faceIdx = 0; faceIdx < faces.size(); ++faceIdx) {
        const Face& face = faces[faceIdx];
        for (int i = 0, j = 2; i < 3; j = i++) {
            if (connectivity.faceFaceTable[faceIdx][j] < 0) {
                // no neighbor -> both vertices boundary
                boundary[face.vertex(i)] = true;
                boundary[face.vertex(j)] = true;
            }
        }
    }
}

void Mesh::computeNormals() {
    normals.resize(faces.size());
    for (uint i = 0; i < faces.size(); ++i) {
        const Face& f = faces[i];
        const math::Point3 dp1 = vertices[f.b] - vertices[f.a];
        const math::Point3 dp2 = vertices[f.c] - vertices[f.a];
        const math::Point3 n = math::crossProduct(dp1, dp2);
        const double length = math::ublas::norm_2(n);
        if (length > 0.) {
            normals[i] = n * (1.0 / length);
        } else {
            // degenerated face, assign arbitrary normal
            normals[i] = { 0, 0, 1 };
        }
    }
}

void Mesh::computeVertexFaceTable() {
    std::vector<std::vector<int>>& vf = connectivity.vertexFaceTable;
    vf.resize(vertices.size());

    // prepare lists: average vertex degree is 6
    for (std::vector<int>& vec : vf) {
        vec.clear();
        vec.reserve(8);
    }

    for (uint i = 0; i < faces.size(); i++) {
        const Face& face = faces[i];
        for (uint j = 0; j < 3; ++j) {
            vf[face.vertex(j)].push_back(i);
        }
    }
}


void Mesh::computeVertexVertexTable() {
    std::vector<std::vector<int>>& vv = connectivity.vertexVertexTable;
    vv.resize(vertices.size());
    for (std::vector<int>& vec : vv) {
        vec.clear();
        vec.reserve(16);
    }

    for (const Face& face : faces) {
        int a = face.a, b = face.b, c = face.c;

        vv[a].push_back(b);
        vv[a].push_back(c);
        vv[b].push_back(c);
        vv[b].push_back(a);
        vv[c].push_back(a);
        vv[c].push_back(b);
    }

    // remove duplicities in lists
    for (std::vector<int>& vec : vv) {
        std::sort(vec.begin(), vec.end());
        auto last = std::unique(vec.begin(), vec.end());
        vec.erase(last, vec.end());
        vec.shrink_to_fit();
    }
}

namespace {
struct EdgeKey {
    int v1, v2; // vertex indices

    EdgeKey(const int a, const int b)
        : v1(a)
        , v2(b) {
        if (v1 > v2) {
            std::swap(v1, v2);
        }
    }

    bool operator<(const EdgeKey& other) const {
        return std::make_tuple(v1, v2) < std::make_tuple(other.v1, other.v2);
    }
};

struct EdgeData {
    int f1, f2; // face indices

    EdgeData()
        : f1(-1)
        , f2(-1) {}

    void addFace(int f) {
        if (f1 < 0) {
            f1 = f;
        } else if (f2 < 0) {
            f2 = f;
        } else {
            LOG(warn1) << "Mesh warning: edge shared by more than two triangles.";
        }
    }

    int getNeighbor(int f) const {
        return (f == f1) ? f2 : f1;
    }
};
}

void Mesh::computeFaceFaceTable() {
    std::map<EdgeKey, EdgeData> edges;

    for (uint i = 0; i < faces.size(); i++) {
        const Face& face = faces[i];
        edges[EdgeKey(face.a, face.b)].addFace(i);
        edges[EdgeKey(face.b, face.c)].addFace(i);
        edges[EdgeKey(face.c, face.a)].addFace(i);
    }

    connectivity.faceFaceTable.resize(faces.size());
    UTILITY_OMP(parallel for)
    for (uint i = 0; i < faces.size(); i++) {
        const Face& face = faces[i];
        std::array<int, 3>& ni = connectivity.faceFaceTable[i];
        ni[0] = edges[EdgeKey(face.a, face.b)].getNeighbor(i);
        ni[1] = edges[EdgeKey(face.b, face.c)].getNeighbor(i);
        ni[2] = edges[EdgeKey(face.c, face.a)].getNeighbor(i);
    }
}

void Mesh::recomputeTopologyData() {
    // order matters here
    if (!connectivity.faceFaceTable.empty()) {
        computeFaceFaceTable();
    }
    if (!boundary.empty()) {
        computeBoundary();
    }
    if (!connectivity.vertexFaceTable.empty()) {
        computeVertexFaceTable();
    }
    if (!connectivity.vertexVertexTable.empty()) {
        computeVertexVertexTable();
    }
}

void Mesh::skirt( const math::Point3 & down ) {
    typedef std::size_t Index;

    enum class Status {
            FW, BW, BI
    };

    struct Edge {
        Index v1, v2;
        Index t1, t2;
        mutable Status status;

        Edge(Index pv1, Index pv2, Index pt1, Index pt2)
            : v1( std::min( pv1, pv2 ) ), v2( std::max( pv1, pv2 ) ),
              t1( pv1 <= pv2 ? pt1 : pt2 ), t2( pv1 <= pv2 ? pt2 : pt1 ),
              status( pv1 <= pv2 ? Status::FW : Status::BW ) {}

        void update( Index pv1, Index pv2 ) const {

            if ( pv1 <= pv2 && status == Status::BW ) status = Status::BI;
            if ( pv1 > pv2 && status == Status::FW ) status = Status::BI;
        }


        bool operator < ( const Edge & edge ) const {

            return v1 < edge.v1 || ( v1 == edge.v1 && v2 < edge.v2 );
        }

    };


    typedef std::set<Edge> Edges;
    typedef std::map<Index, Index> DownMap;

    Edges edges;
    DownMap vdownmap, tdownmap;

    // find odd edges
    for ( Face f : faces ) {

        edges.insert( Edge(f.a,f.b,f.ta,f.tb) ).first->update(f.a,f.b);
        edges.insert( Edge(f.b,f.c,f.tb,f.tc) ).first->update(f.b,f.c);
        edges.insert( Edge(f.c,f.a,f.tc,f.ta) ).first->update(f.c,f.a);
    }

    // iterate through edges
    int evenc(0), oddc(0);

    for ( Edge edge : edges )
        if ( edge.status != Status::BI ) {

            // add new vertexes and tcoords
            if ( vdownmap.find( edge.v1 ) == vdownmap.end() ) {

                vertices.push_back( vertices[edge.v1]  + down );
                vdownmap[edge.v1] = vertices.size()-1;
            }

            if ( tdownmap.find( edge.t1 ) == tdownmap.end() ) {

                tCoords.push_back( tCoords[edge.t1] );
                tdownmap[edge.t1] = tCoords.size()-1;
            }

            if ( vdownmap.find( edge.v2 ) == vdownmap.end() ) {

                vertices.push_back( vertices[edge.v2]  + down );
                vdownmap[edge.v2] = vertices.size()-1;
            }

            if ( tdownmap.find( edge.t2 ) == tdownmap.end() ) {

                tCoords.push_back( tCoords[edge.t2] );
                tdownmap[edge.t2] = tCoords.size()-1;
            }

            // add new faces
            if ( edge.status == Status::FW ) {

                faces.emplace_back( vdownmap[edge.v1], edge.v2, edge.v1,
                           tdownmap[edge.t1], edge.t2, edge.t1 );
                faces.emplace_back( vdownmap[edge.v1], vdownmap[edge.v2], edge.v2,
                           tdownmap[edge.t1], tdownmap[edge.t2], edge.t2 );
            }

            if ( edge.status == Status::BW ) {

                faces.emplace_back( vdownmap[edge.v1], edge.v1, edge.v2,
                           tdownmap[edge.t1], edge.t1, edge.t2 );
                faces.emplace_back( vdownmap[edge.v1], edge.v2, vdownmap[edge.v2],
                           tdownmap[edge.t1], edge.t2, tdownmap[edge.t2] );
            }

            oddc++;

        } else {

            evenc++;
        }

    LOG( info1 ) << evenc << " even, " << oddc << " odd.";
}

double Mesh::area(const Face &face) const
{
    return (norm_2(math::crossProduct(b(face) - a(face)
                                     , c(face) - a(face)))
            * 0.5);
}

double Mesh::txArea(const Face &face) const
{
    return (std::abs(math::crossProduct(math::Point2(tb(face) - ta(face))
                                        , math::Point2(tc(face) - ta(face)))
                     * 0.5));
}

math::Point3 Mesh::barycenter(const Face &face) const
{
    const auto &va(a(face));
    const auto &vb(b(face));
    const auto &vc(c(face));

    return {
        (va(0) + vb(0) + vc(0)) / 3.0
        , (va(1) + vb(1) + vc(1)) / 3.0
        , (va(2) + vb(2) + vc(2)) / 3.0
    };
}

} //namespace geometry

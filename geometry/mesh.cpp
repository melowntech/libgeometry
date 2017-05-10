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

void Mesh::skirt( const math::Point3 & down ) {

    (void) down;

    enum class Status {
            FW, BW, BI
    };

    struct Edge {
        int v1, v2;
        int t1, t2;
        mutable Status status;

        Edge(int pv1, int pv2, int pt1, int pt2)
            : v1( std::min( pv1, pv2 ) ), v2( std::max( pv1, pv2 ) ),
              t1( pv1 <= pv2 ? pt1 : pt2 ), t2( pv1 <= pv2 ? pt2 : pt1 ),
              status( pv1 <= pv2 ? Status::FW : Status::BW ) {}

        void update( int pv1, int pv2 ) const {

            if ( pv1 <= pv2 && status == Status::BW ) status = Status::BI;
            if ( pv1 > pv2 && status == Status::FW ) status = Status::BI;
        }


        bool operator < ( const Edge & edge ) const {

            return v1 < edge.v1 || ( v1 == edge.v1 && v2 < edge.v2 );
        }

    };


    typedef std::set<Edge> Edges;
    typedef std::map<int,int> DownMap;

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

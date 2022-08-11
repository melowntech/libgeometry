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


namespace geometry {

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

                faces.emplace_back(
                    static_cast<geometry::Face::index_type>(vdownmap[edge.v1])
                    , static_cast<geometry::Face::index_type>(edge.v2)
                    , static_cast<geometry::Face::index_type>(edge.v1)
                    , static_cast<geometry::Face::index_type>(tdownmap[edge.t1])
                    , static_cast<geometry::Face::index_type>(edge.t2)
                    , static_cast<geometry::Face::index_type>(edge.t1));
                faces.emplace_back(
                    static_cast<geometry::Face::index_type>(vdownmap[edge.v1])
                    , static_cast<geometry::Face::index_type>(vdownmap[edge.v2])
                    , static_cast<geometry::Face::index_type>(edge.v2)
                    , static_cast<geometry::Face::index_type>(tdownmap[edge.t1])
                    , static_cast<geometry::Face::index_type>(tdownmap[edge.t2])
                    , static_cast<geometry::Face::index_type>(edge.t2));
            }

            if ( edge.status == Status::BW ) {

                faces.emplace_back(
                    static_cast<geometry::Face::index_type>(vdownmap[edge.v1])
                    , static_cast<geometry::Face::index_type>(edge.v1)
                    , static_cast<geometry::Face::index_type>(edge.v2)
                    , static_cast<geometry::Face::index_type>(tdownmap[edge.t1])
                    , static_cast<geometry::Face::index_type>(edge.t1)
                    , static_cast<geometry::Face::index_type>(edge.t2));
                faces.emplace_back(
                    static_cast<geometry::Face::index_type>(vdownmap[edge.v1])
                    , static_cast<geometry::Face::index_type>(edge.v2)
                    , static_cast<geometry::Face::index_type>(vdownmap[edge.v2])
                    , static_cast<geometry::Face::index_type>(tdownmap[edge.t1])
                    , static_cast<geometry::Face::index_type>(edge.t2)
                    , static_cast<geometry::Face::index_type>(
                        tdownmap[edge.t2]));
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

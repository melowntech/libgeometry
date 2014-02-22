#include <stdexcept>
#include <fstream>

#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>

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

Mesh::pointer Mesh::simplify( int faceCount ) const
{
    OMMesh omMesh;
    toOpenMesh(*this, omMesh);

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



void Mesh::convertTo( geometry::Obj & obj ) const {
    
    for ( math::Point3 vertex : vertices ) {
        
        geometry::Obj::Vector3d overtex;
        overtex.x = vertex[0]; overtex.y = vertex[1]; overtex.z = vertex[2];
        
        obj.addVertex( overtex );
        
        //LOG( info1 ) << "[" << overtex.x << ", " << overtex.y << ", " 
        //    << overtex.z << "]";
    }

    for ( math::Point2 texture : tCoords ) {
        
        geometry::Obj::Vector3d otexture;
        otexture.x = texture[0]; otexture.y = texture[1]; otexture.z = 0.0;
        
        obj.addTexture( otexture );
    }
    
    for ( geometry::Face face : faces ) {
        
        geometry::Obj::Facet facet;
        
        facet.v[0] = face.a; facet.v[1] = face.b; facet.v[2] = face.c;
        facet.t[0] = face.ta; facet.t[1] = face.tb; facet.t[2] = face.tc;
        
        obj.addFacet( facet );
    }
}


} //namespace geometry

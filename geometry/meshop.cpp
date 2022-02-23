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
/**
 * @file math.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * 3D mesh operations
 */

#include <unordered_set>
#include <unordered_map>

#include "meshop.hpp"
#include "parse-obj.hpp"
#include "parse-ply.hpp"
#include "triclip.hpp"

#include "utility/expect.hpp"
#include "utility/small_list.hpp"

#include <set>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/lexical_cast.hpp>

namespace geometry {

namespace ublas = boost::numeric::ublas;
namespace ba = boost::algorithm;

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

        facet.v[0] = int(face.a);
        facet.v[1] = int(face.b);
        facet.v[2] = int(face.c);

        facet.t[0] = int(face.ta);
        facet.t[1] = int(face.tb);
        facet.t[2] = int(face.tc);

        obj.addFacet( facet );
    }

    return obj;
}

Mesh::pointer asMesh(const Obj &obj){
    auto newMesh(std::make_shared<geometry::Mesh>());
    for( const auto&v : obj.vertices ){
        newMesh->vertices.push_back(v);
    }

    for( const auto&t : obj.texcoords ){
        newMesh->tCoords.emplace_back( t(0), t(1) );
    }

    for( const auto&f : obj.facets ){
        newMesh->addFace(f.v[0], f.v[1], f.v[2], f.t[0], f.t[1], f.t[2]);
    }

    return newMesh;
}

std::string ObjMaterial::name(std::size_t index) const
{
    if (index < names.size()) { return names[index]; }
    return boost::lexical_cast<std::string>(index);
}

namespace {

void addMtl(std::ostream &out, const ObjMaterial &mtl, unsigned int imageId)
{
    if (mtl.libs.empty()) { return; }
    out << "usemtl " << mtl.name(imageId) << '\n';
}

} // namespace

void saveAsObj(const Mesh &mesh, std::ostream &out
               , const ObjMaterial &mtl
               , const boost::filesystem::path &filepath
               , const ObjStreamSetup &streamSetup)
{
    // keep current numeric precision for future use
    const auto oldPrecision(out.precision());

    for (const auto &lib : mtl.libs) {
        out << "mtllib " << lib << '\n';
    }

    if (!streamSetup.vertex(out)) {
        out.setf(std::ios::scientific, std::ios::floatfield);
    }

    for (const auto &vertex : mesh.vertices) {
        out << "v " << vertex(0) << ' ' << vertex(1) << ' '  << vertex(2)
            << '\n';
    }

    if (!streamSetup.tx(out)) {
        out.setf(std::ios::scientific, std::ios::floatfield);
        // reset precision to recorded one
        out.precision(oldPrecision);
    }

    for (const auto &tCoord : mesh.tCoords) {
        out << "vt " << tCoord(0) << ' ' << tCoord(1) << '\n';
    }

    const bool textured(!mesh.tCoords.empty());

    unsigned int currentImageId(static_cast<unsigned int>(-1));

    for (const auto &face : mesh.faces) {
        if (face.degenerate()) {
            continue;
        }

        if (textured && (face.imageId != currentImageId)) {
            addMtl(out, mtl, face.imageId);
            currentImageId = face.imageId;
        }

        if (textured) {
            out << "f " << face.a + 1 << '/' << face.ta + 1 << "/ "
                << face.b + 1 << '/' << face.tb + 1 << "/ "
                << face.c + 1 << '/' << face.tc + 1 << "/\n";
        } else {
            out << "f " << face.a + 1 << ' '
                << face.b + 1 << ' '
                << face.c + 1 << "\n";
        }
    }

    if (!out) {
        LOGTHROW(err3, std::runtime_error)
            << "Unable to save mesh to <" << filepath << ">.";
    }
}

void saveAsObj(const Mesh &mesh, std::ostream &out
               , const ObjMaterial &mtl
               , const boost::filesystem::path &filepath
               , bool setFormat)
{
    struct DontSetFormat : ObjStreamSetup {
        virtual bool vertex(std::ostream&) const { return true; }
        virtual bool tx(std::ostream&) const { return true; }
    };

    const ObjStreamSetup &streamSetup(setFormat
                                      ? ObjStreamSetup()
                                      : DontSetFormat());

    saveAsObj(mesh, out, mtl, filepath, streamSetup);
}

void saveAsObj(const Mesh &mesh, const boost::filesystem::path &filepath
               , const ObjMaterial &mtl, const ObjStreamSetup &streamSetup)
{
    LOG(info2) << "Saving mesh to file <" << filepath << ">.";

    std::ofstream f;
    f.exceptions(std::ios::badbit | std::ios::failbit);
    try {
        f.open(filepath.string(), std::ios_base::out | std::ios_base::trunc);
    } catch (const std::exception&) {
        LOGTHROW(err3, std::runtime_error)
            << "Unable to save mesh to <" << filepath << ">.";
    }

    saveAsObj(mesh, f, mtl, filepath, streamSetup);
}

void saveAsPly( const Mesh &mesh, const boost::filesystem::path &filepath){
    LOG(info2) << "Saving mesh to file <" << filepath << ">.";

    std::ofstream out;
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out.open(filepath.string(), std::ios_base::out | std::ios_base::trunc);
    out.setf(std::ios::scientific, std::ios::floatfield);

    unsigned validFaces(0);
    for (const auto &face : mesh.faces) {
        if (!face.degenerate() && mesh.good(face))
            validFaces++;
    }

    out << "ply\n"
        << "format ascii 1.0\n"
        << "comment generated by window-mesh\n"
        << "element vertex " << mesh.vertices.size() << '\n'
        << "property float x\n"
        << "property float y\n"
        << "property float z\n"
        << "element face " << validFaces << '\n'
        << "property list uchar int vertex_indices\n"
        << "end_header\n";

    for (const auto &vertex : mesh.vertices)
    {
        out << vertex(0) << ' ' << vertex(1) << ' '  << vertex(2) << '\n';
    }

    for (const auto &face : mesh.faces)
    {
        if (face.degenerate()) {
            continue;
        }
        if (!mesh.good(face)) {
            LOG(warn2) << "Invalid vertex index in face.";
            continue;
        }
        out << "3 " << face.a << ' ' << face.b << ' ' << face.c << '\n';
    }

    out.close();

    if (!out) {
        LOGTHROW(err3, std::runtime_error)
            << "Unable to save mesh to <" << filepath << ">.";
    }
}

Mesh loadPly(const boost::filesystem::path& filepath)
{
    class SimplePlyParser : public PlyParserBase
    {
    public:
        Mesh mesh;

        void loadHeader(const std::vector<PlyElement>& elements,
                        const std::vector<std::string>& /* comments */) override
        {
            for (const auto& el : elements)
            {
                if (el.name == "vertex")
                {
                    mesh.vertices.reserve(el.count);

                    // check vertex properties
                    if (el.props.size() != 3)
                    {
                        LOGTHROW(err4, std::runtime_error)
                            << "No additional vertex properties supported in "
                               "simple parser";
                    }
                    utility::expect(el.props[0].name == "x");
                    utility::expect(el.props[1].name == "y");
                    utility::expect(el.props[2].name == "z");
                }
                else if (el.name == "face")
                {
                    mesh.faces.reserve(el.count);

                    // check face properties
                    if (el.props.size() != 1)
                    {
                        LOGTHROW(err4, std::runtime_error)
                            << "No additional face properties supported in "
                               "simple parser";
                    }
                    if (!ba::starts_with(el.props[0].type, "list"))
                    {
                        LOGTHROW(err4, std::runtime_error)
                            << "Expected face property type to be a list, but "
                               "got: "
                            << el.props[0].type;
                    }
                    // does not check the name which should generally be
                    // "vertex_index" but might vary
                }
                else
                {
                    LOGTHROW(err4, std::runtime_error)
                        << "Unexpected element in header - not supported by "
                           "simple parser: "
                        << el.name;
                }
            }
        }

        void addVertex(std::istream& f) override
        {
            double x, y, z;
            f >> x >> y >> z;
            mesh.vertices.emplace_back(x, y, z);
        }

        void addFace(std::istream& f) override
        {
            int n, a, b, c;
            f >> n >> a >> b >> c;
            if (n != 3)
            {
                LOGTHROW(err4, std::runtime_error)
                    << "Simple parser only supports loading triangular meshes, "
                       "but got face with " << n << " vertices.";
            }
            mesh.faces.emplace_back(a, b, c);
        }
    } simpleParser;

    parsePly(simpleParser, filepath);

    LOG(info1) << "Loaded mesh with " << simpleParser.mesh.vertices.size()
               << " vertices and " << simpleParser.mesh.faces.size()
               << " faces";

    return simpleParser.mesh;
}

const unsigned int imageIdLimit(1<<16);

void loadObj( ObjParserBase &parser
            , const boost::filesystem::path &filename)
{
    std::ifstream file;
    file.exceptions(std::ios::badbit | std::ios::failbit);
    file.open(filename.string());

    if (!parser.parse(file)) {
        LOGTHROW(err2, std::runtime_error)
            << "LoadObj failed to parse file " << filename << ".";
    }
}

Mesh loadObj(const boost::filesystem::path &filename, ObjMaterial *mtl)
{
    struct Obj2MeshParser : public ObjParserBase {
        Obj2MeshParser() : imageId(), namedCount() {}

        void addVertex( const Vector3d &v ) {
            mesh.vertices.push_back(v);
        }

        void addTexture( const Vector3d &t ) {
            mesh.tCoords.emplace_back(t.x, t.y);
        }

        void addFacet( const Facet &f ) {
            mesh.addFace( f.v[0], f.v[1], f.v[2]
                        , f.t[0], f.t[1], f.t[2] );
            mesh.faces.back().imageId = imageId;
        }

        void addNormal( const Vector3d& ) { }

        void materialLibrary(const std::string &lib) { seenLibs.insert(lib); }
        void useMaterial(const std::string &name) {
            auto fmtlMap(mtlMap.find(name));
            if (fmtlMap == mtlMap.end()) {
                unsigned int id(0);
                try {
                    // is id a (non negative) number?
                    auto number(boost::lexical_cast<unsigned int>(name));
                    id = number;
                } catch (const boost::bad_lexical_cast&) {
                    // not a number
                    id = imageIdLimit + namedCount++;
                }

                fmtlMap = mtlMap.insert(MtlMap::value_type(name, id)).first;
            }

            imageId = fmtlMap->second;
        }

        Mesh mesh;
        unsigned int imageId;
        std::set<std::string> seenLibs;
        typedef std::map<std::string, int> MtlMap;
        MtlMap mtlMap;
        int namedCount;
    } parser_;

    loadObj(parser_, filename);

    std::vector<std::string> materials;

    // regenerate IDs if there was any non-numeric material
    if (parser_.namedCount) {
        typedef std::map<unsigned int, unsigned int> Mapping;
        Mapping mapping;
        for (const auto &item : parser_.mtlMap) {
            mapping[item.second] = materials.size();
            materials.push_back(item.first);
        }

        for (auto &face : parser_.mesh.faces) {
            face.imageId = mapping[face.imageId];
        }
    }

    if (mtl) {
        mtl->libs.assign(parser_.seenLibs.begin(), parser_.seenLibs.end());
        mtl->names = std::move(materials);
    }

    return parser_.mesh;
}

namespace {

void clipImpl(const Mesh &omesh, Mesh &mesh
              , const std::vector<ClipPlane> &planes)
{
    ClipTriangle::list clipped;
    for (const auto &face : omesh.faces) {
        clipped.emplace_back(
              omesh.vertices[face.a]
            , omesh.vertices[face.b]
            , omesh.vertices[face.c]);
    }

    std::vector<double> tinfos;
    for (const auto &plane : planes) {
        clipped = clipTriangles(clipped, plane, tinfos);
    }

    typedef math::Points3::size_type Index;
    std::map<math::Point3, Index> pMap;
    Index next = 0;

    for (const auto &triangle : clipped)
    {
        Index indices[3];
        for (int i = 0; i < 3; i++)
        {
            auto pair = pMap.insert(std::make_pair(triangle.pos[i], next));
            if (pair.second) {
                next++;
            }
            indices[i] = pair.first->second;

            if (indices[i] >= mesh.vertices.size()) {
                mesh.vertices.push_back(triangle.pos[i]);
            }
        }
        // do not add degenerated faces
        if ((indices[0] != indices[1]) &&
            (indices[1] != indices[2]) &&
            (indices[0] != indices[2]))
        {
            mesh.addFace(indices[0], indices[1], indices[2]);
        }
    }
}

std::vector<ClipPlane> planes(const math::Extents2 &extents)
{
    std::vector<ClipPlane> planes;
    planes.emplace_back(+1.,  0., 0., extents.ll[0]);
    planes.emplace_back(-1.,  0., 0., -extents.ur[0]);
    planes.emplace_back(0.,  +1., 0., extents.ll[1]);
    planes.emplace_back(0.,  -1., 0., -extents.ur[1]);
    return planes;
}

std::vector<ClipPlane> planes(const math::Extents3 &extents)
{
    std::vector<ClipPlane> planes;
    planes.emplace_back(+1.,  0., 0., extents.ll[0]);
    planes.emplace_back(-1.,  0., 0., -extents.ur[0]);
    planes.emplace_back(0.,  +1., 0., extents.ll[1]);
    planes.emplace_back(0.,  -1., 0., -extents.ur[1]);
    planes.emplace_back(0.,  0., +1., extents.ll[2]);
    planes.emplace_back(0.,  0., -1., -extents.ur[2]);
    return planes;
}

} // namespace

Mesh clip(const Mesh &omesh, const math::Extents2 &extents)
{
    Mesh out;
    clipImpl(omesh, out, planes(extents));
    return out;
}

Mesh clip(const Mesh &omesh, const math::Extents3 &extents)
{
    Mesh out;
    clipImpl(omesh, out, planes(extents));
    return out;
}

Mesh::pointer removeNonManifoldEdges(Mesh omesh)
{
    auto ofaces = omesh.faces;
    auto pmesh(std::make_shared<geometry::Mesh>(std::move(omesh)));
    auto& mesh(*pmesh);
    mesh.faces.clear();
    mesh.faces.shrink_to_fit();

    typedef Face::index_type index_type;

    struct EdgeKey
    {
        index_type v1, v2; // vertex indices

        EdgeKey(index_type v1, index_type v2)
        {
            this->v1 = std::min(v1, v2);
            this->v2 = std::max(v1, v2);
        }

        bool operator< (const EdgeKey& other) const
        {
            return (v1 == other.v1) ? (v2 < other.v2) : (v1 < other.v1);
        }
    };

    struct Edge {
        utility::small_list<index_type, 2> facesIndices;
    };

    //count faces for each edge
    std::map<EdgeKey,Edge> edgeMap;
    for(index_type fi=0; fi<ofaces.size(); fi++){
        const auto & face(ofaces[fi]);
        EdgeKey edgeKeys[3] = { EdgeKey(face.a,face.b)
                           , EdgeKey(face.b,face.c)
                           , EdgeKey(face.c,face.a) };
        for(const auto & key : edgeKeys){
            auto it = edgeMap.find(key);
            if(it==edgeMap.end()){
                //if edge is not present insert it with current face
                Edge edge;
                edge.facesIndices.insert(fi);
                edgeMap.insert(std::make_pair(key,edge));
            }
            else{
                //if edge is present add current face to it
                it->second.facesIndices.insert(fi);
            }
        }
    }

    //collect faces incident with non-manifold edge
    std::set<index_type> facesToOmit;
    for(auto it = edgeMap.begin(); it!=edgeMap.end(); it++){
        if(it->second.facesIndices.size()>2){
            it->second.facesIndices.for_each([&](index_type fi) {
                facesToOmit.insert(fi);
            });
        }
    }

    for(index_type fi=0; fi<ofaces.size(); fi++){
        const auto & face(ofaces[fi]);
        if(facesToOmit.find(fi)==facesToOmit.end()){
            mesh.addFace( face.a, face.b, face.c
                        , face.ta, face.tb, face.tc );
        }
    }

    return pmesh;
}

Mesh::pointer removeIsolatedVertices( const Mesh& imesh ){
    auto pmesh(std::make_shared<geometry::Mesh>());
    auto & mesh(*pmesh);

    std::map<math::Points3::size_type, math::Points3::size_type> vertexMap;
    std::map<math::Points2::size_type, math::Points2::size_type> tCoordsMap;

    for( const auto& face: imesh.faces ){
        math::Points3::size_type vindices[3] { face.a , face.b, face.c };
        math::Points2::size_type tindices[3] { face.ta , face.tb, face.tc };

        for(unsigned int i=0; i<3; ++i){
            auto vit = vertexMap.find(vindices[i]);
            auto tit = tCoordsMap.find(tindices[i]);

            if(vit == vertexMap.end()){
                mesh.vertices.push_back(imesh.vertices[vindices[i]]);
                vit = vertexMap.insert(std::make_pair(vindices[i],mesh.vertices.size()-1)).first;
            }
            if(imesh.tCoords.size() > 0 && tit == tCoordsMap.end()){
                mesh.tCoords.push_back(imesh.vertices[tindices[i]]);
                tit = tCoordsMap.insert(std::make_pair(tindices[i],mesh.tCoords.size()-1)).first;
            }

            vindices[i] = vit->second;
            tindices[i] = tit->second;
        }

        if(imesh.tCoords.size() > 0){
            mesh.addFace( vindices[0], vindices[1], vindices[2]
                        , tindices[0], tindices[1], tindices[2] );
        }
        else{
            mesh.addFace( vindices[0], vindices[1], vindices[2] );
        }
    }

    return pmesh;
}


Mesh::pointer refine( const Mesh & omesh, unsigned int maxFacesCount)
{
    auto pmesh(std::make_shared<geometry::Mesh>(omesh));
    auto & mesh(*pmesh);

    struct EdgeKey
    {
        std::size_t v1, v2; // vertex indices

        EdgeKey(std::size_t v1, std::size_t v2)
        {
            this->v1 = std::min(v1, v2);
            this->v2 = std::max(v1, v2);
        }

        bool operator< (const EdgeKey& other) const
        {
            return (v1 == other.v1) ? (v2 < other.v2) : (v1 < other.v1);
        }
    };

    struct Edge {
        typedef enum {
            AB,
            BC,
            CA
        } EdgeType;

        std::size_t v1, v2;
        int f1, f2;
        EdgeType et1, et2;

        float length;

        Edge(std::size_t pv1, std::size_t pv2, float plength)
            :v1(std::min(pv1,pv2)),v2(std::max(pv1,pv2)), length(plength){
            f1 = -1;
            f2 = -1;
        }

        void addFace(std::size_t pv1, std::size_t pv2, int fid, EdgeType type)
        {
            if(pv1<pv2){
                f1=fid;
                et1 = type;
            }
            else{
                f2=fid;
                et2 = type;
            }
        }

        bool operator< (const Edge& other) const
        {
            return length < other.length;
        }
    };

    struct EdgeMap {
        std::map<EdgeKey, std::shared_ptr<Edge>> map;
        std::vector<std::shared_ptr<Edge>> heap;

        bool compareEdgePtr(const std::shared_ptr<Edge> &a, const std::shared_ptr<Edge> &b){
            return a->length<b->length;
        }

        void addFaceEdge(std::size_t pv1, std::size_t pv2, int fid
                         , Edge::EdgeType type, float length)
        {
            EdgeKey key(pv1, pv2);
            auto it(map.find(key));
            if(it!=map.end()){
                it->second->addFace(pv1, pv2, fid, type );
            }
            else{
                heap.push_back(std::make_shared<Edge>(pv1, pv2, length));
                heap.back()->addFace(pv1, pv2, fid, type );
                map.insert(std::make_pair(key,heap.back()));
                std::push_heap(heap.begin(),heap.end(), [this]( const std::shared_ptr<Edge> &a
                                                              , const std::shared_ptr<Edge> &b){
                    return this->compareEdgePtr(a,b);
                });
            }
        }

        Edge pop_top_edge(){
            auto edge = *heap[0];
            std::pop_heap(heap.begin(), heap.end(), [this]( const std::shared_ptr<Edge> &a
                                                          , const std::shared_ptr<Edge> &b){
                return this->compareEdgePtr(a,b);
            });
            heap.pop_back();
            map.erase(EdgeKey(edge.v1, edge.v2));
            return edge;
        }

        Edge top_edge(){
            return *heap[0];
        }

        void addFaceEdges(const Mesh & mesh, int fid){
            const auto& f = mesh.faces[fid];
            auto e1Length =  float(ublas::norm_2(mesh.vertices[f.a]-mesh.vertices[f.b]));
            addFaceEdge(f.a, f.b, fid, Edge::EdgeType::AB, e1Length);

            auto e2Length =  float(ublas::norm_2(mesh.vertices[f.b]-mesh.vertices[f.c]));
            addFaceEdge(f.b, f.c, fid, Edge::EdgeType::BC, e2Length);

            auto e3Length =  float(ublas::norm_2(mesh.vertices[f.c]-mesh.vertices[f.a]));
            addFaceEdge(f.c, f.a, fid, Edge::EdgeType::CA, e3Length);
        }

        std::size_t size(){
            return heap.size();
        }

    };

    EdgeMap edgeMap;

    auto splitEdge = [&mesh,&edgeMap]( int fid, Edge::EdgeType type
                       , std::size_t vid) mutable -> void{
        auto & face = mesh.faces[fid];
        switch(type){
            case Edge::EdgeType::AB:
                {
                    if(mesh.tCoords.size()>0){
                        math::Point2 tcMiddle = ( mesh.tCoords[face.ta]
                                                + mesh.tCoords[face.tb]) * 0.5;
                        mesh.tCoords.push_back(tcMiddle);
                    }
                    mesh.addFace( mesh.faces[fid].b, mesh.faces[fid].c
                                , vid
                                , mesh.faces[fid].tb, mesh.faces[fid].tc
                                , mesh.tCoords.size()-1);

                    mesh.faces[fid].b = vid;
                    mesh.faces[fid].tb = mesh.tCoords.size()-1;

                    edgeMap.addFaceEdges(mesh, fid);
                    edgeMap.addFaceEdges(mesh, int(mesh.faces.size()-1));
                }
                break;
            case Edge::EdgeType::BC:
                {
                    if(mesh.tCoords.size()>0){
                        math::Point2 tcMiddle = (mesh.tCoords[face.tb]
                                                + mesh.tCoords[face.tc]) * 0.5;
                        mesh.tCoords.push_back(tcMiddle);
                    }
                    mesh.addFace( mesh.faces[fid].c,mesh.faces[fid].a, vid
                                , mesh.faces[fid].tc,mesh.faces[fid].ta
                                , mesh.tCoords.size()-1);

                    mesh.faces[fid].c = vid;
                    mesh.faces[fid].tc = mesh.tCoords.size()-1;

                    edgeMap.addFaceEdges(mesh, fid);
                    edgeMap.addFaceEdges(mesh, int(mesh.faces.size()-1));
                }
                break;
            case Edge::EdgeType::CA:
                {
                    if(mesh.tCoords.size()>0){
                        math::Point2 tcMiddle = (mesh.tCoords[face.tc]
                                                + mesh.tCoords[face.ta]) * 0.5;
                        mesh.tCoords.push_back(tcMiddle);
                    }

                    mesh.addFace( mesh.faces[fid].a,mesh.faces[fid].b, vid
                                , mesh.faces[fid].ta,mesh.faces[fid].tb
                                , mesh.tCoords.size()-1);

                    mesh.faces[fid].a = vid;
                    mesh.faces[fid].ta = mesh.tCoords.size()-1;

                    edgeMap.addFaceEdges(mesh, fid);
                    edgeMap.addFaceEdges(mesh, int(mesh.faces.size()-1));
                }
                break;
        }
    };

    for (std::size_t i=0; i<mesh.faces.size(); ++i ) {
        //add all 3 edges
        edgeMap.addFaceEdges(mesh, int(i));
    }

    //sort edges by length
    while( mesh.faces.size() < maxFacesCount && edgeMap.size()>0 ){
        //split edge
        auto edge = edgeMap.pop_top_edge();

        //find middle
        math::Point3 middle = (mesh.vertices[edge.v1]
                            + mesh.vertices[edge.v2]) * 0.5;
        mesh.vertices.push_back(middle);

        //split first face
        if(edge.f1>=0){
            splitEdge(edge.f1, edge.et1, mesh.vertices.size()-1);
        }

        //split second face
        if(edge.f2>=0){
            splitEdge(edge.f2, edge.et2, mesh.vertices.size()-1);
        }
    }

    return pmesh;
}

MeshInfo measurePly(std::istream &is, const boost::filesystem::path &path)
{
    // read header
    std::string line;
    int nvert = -1, ntris = -1;
    do {
        if (getline(is, line).eof()) break;
        std::sscanf(line.c_str(), "element vertex %d", &nvert);
        std::sscanf(line.c_str(), "element face %d", &ntris);
    } while (line != "end_header");

    if (nvert < 0 || ntris < 0) {
        LOGTHROW(err2, std::runtime_error)
            << "Unknown PLY format in file " << path << ".";
    }

    return MeshInfo(nvert, ntris);
}

void loadPly(ObjParserBase &parser, std::istream &is
             , const boost::filesystem::path &path)
{
    const auto mi(measurePly(is, path));

    // load points
    ObjParserBase::Vector3d v;
    for (std::size_t i = 0; i < mi.vertexCount; i++) {
        is >> v.x >> v.y >> v.z;
        parser.addVertex(v);
    }

    // load triangles
    ObjParserBase::Facet f;
    int n;
    for (std::size_t i = 0; i < mi.faceCount; i++) {
        is >> n >> f.v[0] >> f.v[1] >> f.v[2];
        utility::expect(n == 3
                        , "Only triangles are supported in PLY files (&s)."
                        , path);
        parser.addFacet(f);
    }
}

MeshInfo measurePly(const boost::filesystem::path &path)
{
    std::ifstream f(path.string());
    if (!f.good()) {
        LOGTHROW(err2, std::runtime_error)
            << "Can't open PLY file " << path << ".";
    }

    f.exceptions(std::ios::badbit | std::ios::failbit);
    const auto mi(measurePly(f, path));
    f.close();

    return mi;
}

bool objHasVertexOrFace(const boost::filesystem::path &path)
{
    bool hasVorF(false);

    std::ifstream f(path.string());
    if (!f.good()) {
        LOGTHROW(err2, std::runtime_error)
            << "Can't open OBJ file " << path << ".";
    }

    f.exceptions(std::ios::badbit);

    std::string line;
    do {
        getline(f, line);
        if (f.eof()) break;
        if (f.fail()) {
            LOGTHROW(err2, std::runtime_error)
                << "Failbit when reading <" << path << ">";
        }

        if (line.empty()) continue;
        if (line[0] == 'v' || line[0] == 'f') {
            hasVorF = true;
            break;
        }
    } while (true);
    f.close();

    return hasVorF;
}

void loadPly(ObjParserBase &parser, const boost::filesystem::path &path)
{
    std::ifstream f(path.string());
    if (!f.good()) {
        LOGTHROW(err2, std::runtime_error)
            << "Can't open PLY file " << path << ".";
    }

    f.exceptions(std::ios::badbit | std::ios::failbit);

    loadPly(parser, f, path);

    f.close();
}

void append(Mesh &mesh, const Mesh &added)
{
    // remember vertex indexing shift and append vertices
    const auto vShift(mesh.vertices.size());
    mesh.vertices.insert(mesh.vertices.end(), added.vertices.begin()
                         , added.vertices.end());

    // remember texturing coord indexing shift and append texturing coordinates
    const auto tcShift(mesh.tCoords.size());
    mesh.tCoords.insert(mesh.tCoords.end(), added.tCoords.begin()
                        , added.tCoords.end());

    // append faces with shifts applied
    for (const auto &f : added.faces) {
        mesh.addFace(f.a + vShift, f.b + vShift, f.c + vShift
                     , f.ta + tcShift, f.tb + tcShift, f.tc + tcShift
                     , f.imageId);
    }
}

Mesh::list splitById(const Mesh &mesh)
{
    typedef unsigned int ImageId;

    struct MeshBuilder {
        typedef std::unordered_map<ImageId, MeshBuilder> map;

        typedef Face::index_type Index;

        typedef std::unordered_map<Index, Index> VMap;

        VMap vertices;
        Mesh mesh;

        Index vertex(const Mesh &m, Index v) {
            auto fvertices(vertices.find(v));
            if (fvertices == vertices.end()) {
                auto nv(mesh.vertices.size());
                mesh.vertices.push_back(m.vertices[v]);
                fvertices = vertices.insert(VMap::value_type(v, nv)).first;
            }
            return fvertices->second;
        }

        void add(const Mesh &m, const Face &face) {
            mesh.faces.emplace_back(vertex(m, face.a)
                                    , vertex(m, face.b)
                                    , vertex(m, face.c)
                                    , 0, 0, 0, face.imageId);
        }
    };

    MeshBuilder::map builders;
    std::unordered_set<ImageId> unique;
    for (const auto &face : mesh.faces) {
        auto fbuilders(builders.find(face.imageId));
        if (fbuilders == builders.end()) {
            fbuilders
                = builders.insert(MeshBuilder::map::value_type
                                  (face.imageId, MeshBuilder())).first;
        }
        fbuilders->second.add(mesh, face);
    }

    Mesh::list out;
    out.reserve(builders.size());
    for (auto &builder : builders) {
        out.push_back(std::move(builder.second.mesh));
    }
    return out;
}

} // namespace geometry

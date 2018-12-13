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

#include <numeric>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <OpenMesh/Tools/Decimater/ModNormalFlippingT.hh>
#include <OpenMesh/Tools/Decimater/ModAspectRatioT.hh>
#include <OpenMesh/Tools/Decimater/ModEdgeLengthT.hh>
#include <boost/numeric/ublas/vector.hpp>

#include "math/math.hpp"

#include "./meshop.hpp"
#include "./detail/hybrid-decimater.hpp"

namespace geometry {

/* code moved from the window-mesh utility */
namespace ublas = boost::numeric::ublas;

namespace {

struct NormalTraits : public OpenMesh::DefaultTraits
{
  FaceAttributes( OpenMesh::Attributes::Normal );
};

typedef OpenMesh::TriMesh_ArrayKernelT<NormalTraits> OMMesh;
typedef OpenMesh::Decimater::DecimaterT<OMMesh> Decimator;
typedef OpenMesh::Decimater::ModQuadricT<OMMesh>::Handle HModQuadric;
typedef detail::ModQuadricConvexT<OMMesh>::Handle HModQuadricConvex;
typedef detail::ModQuadricHybrid<OMMesh>::Handle HModQuadricHybrid;
typedef OpenMesh::Decimater::ModNormalFlippingT<OMMesh>::Handle HModNormalFlippingT;
typedef OpenMesh::Decimater::ModAspectRatioT<OMMesh>::Handle HModAspectRatioT;
typedef OpenMesh::Decimater::ModEdgeLengthT<OMMesh>::Handle HModEdgeLengthT;

/** Convert mesh to OpenMesh data structure and return mesh extents (2D bounding
 *  box)
 */
math::Extents2 toOpenMesh(const geometry::Mesh &mesh, OMMesh& omMesh)
{
    math::Extents2 e(math::InvalidExtents{});

    // create OpenMesh vertices
    std::vector<OMMesh::VertexHandle> handles;
    handles.reserve(mesh.vertices.size());
    for (const auto& v : mesh.vertices) {
        handles.emplace_back(
            omMesh.add_vertex(OMMesh::Point(v(0), v(1), v(2))) );
        update(e, v);
    }

    // create OpenMesh faces
    for (const geometry::Face& face : mesh.faces)
    {
        OMMesh::VertexHandle face_handles[3] = {
            handles[face.a], handles[face.b], handles[face.c]
        };
        omMesh.add_face(face_handles, 3);
    }

    return e;
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

void lockBorder(OMMesh &omMesh, bool inner, bool outer){
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

    double eps = std::max(box[0][1]-box[0][0], box[1][1]-box[1][0])/(1<<10);
    double epsNormal = 0.001;

    for (auto e_it = omMesh.edges_begin();
              e_it != omMesh.edges_end();  ++e_it)
    {
        if(omMesh.is_boundary(e_it))
        {
            //LOG(info3)<<"--------";
            auto hfhandle(omMesh.halfedge_handle(e_it,0));

            auto from(omMesh.from_vertex_handle(hfhandle));
            auto to(omMesh.to_vertex_handle(hfhandle));

            //check if both points lie on the same border
            auto p1(omMesh.point(from));
            auto p2(omMesh.point(to));
            math::Point3 curNormal(
                math::normalize(math::Point3(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])));
            //check if previous edge or next is in the same direction
            bool dirChangeFrom(true);
            bool dirChangeTo(true);
            for(auto vih_it = omMesh.vih_iter(to); vih_it; ++vih_it) {
                if(omMesh.is_boundary(vih_it)){
                    //auto prevhf = omMesh.prev_halfedge_handle(hfhandle);
                    auto pfrom(omMesh.from_vertex_handle(vih_it));
                    auto pto(omMesh.to_vertex_handle(vih_it));
                    auto pp1(omMesh.point(pfrom));
                    auto pp2(omMesh.point(pto));
                    math::Point3 prevNormal(
                    math::normalize(math::Point3(pp2[0]-pp1[0], pp2[1]-pp1[1], pp2[2]-pp1[2])));
                    dirChangeTo = std::abs(ublas::inner_prod(prevNormal,curNormal))<(1.0f-epsNormal);
                }
            }
            for(auto voh_it = omMesh.voh_iter(from); voh_it; ++voh_it) {
                if(omMesh.is_boundary(voh_it)){
                    //auto nexthf = omMesh.next_halfedge_handle(hfhandle);
                    auto nfrom(omMesh.from_vertex_handle(voh_it));
                    auto nto(omMesh.to_vertex_handle(voh_it));
                    auto np1(omMesh.point(nfrom));
                    auto np2(omMesh.point(nto));
                    math::Point3 nextNormal(
                        math::normalize(math::Point3(np2[0]-np1[0], np2[1]-np1[1], np2[2]-np1[2])));
                    dirChangeFrom = std::abs(ublas::inner_prod(nextNormal,curNormal))<(1.0f-epsNormal);
                }
            }

            bool p1Border[4];
            bool p2Border[4];

            p1Border[0] = std::abs(p1[0]-box[0][0])<eps;
            p1Border[1] = std::abs(p1[0]-box[0][1])<eps;
            p1Border[2] = std::abs(p1[1]-box[1][0])<eps;
            p1Border[3] = std::abs(p1[1]-box[1][0])<eps;

            p2Border[0] = std::abs(p2[0]-box[0][0])<eps;
            p2Border[1] = std::abs(p2[0]-box[0][1])<eps;
            p2Border[2] = std::abs(p2[1]-box[1][0])<eps;
            p2Border[3] = std::abs(p2[1]-box[1][0])<eps;

            bool onBorder = false;
            for (int i = 0; i < 4; i++) {
                onBorder = (onBorder || (p1Border[i]==p2Border[i] && p1Border[i]));
            }

            if((inner && !onBorder) || (outer && onBorder)){
                //if they don't lie on the border, lock them
                if(dirChangeFrom){
                    omMesh.status(from).set_locked(true);
                }

                if(dirChangeTo){
                    omMesh.status(to).set_locked(true);
                }
            }
        }
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
            box[i][0] = std::min(double(pt[i]), box[i][0]);
            box[i][1] = std::max(double(pt[i]), box[i][1]);
        }
    }

    // calculate smallest distances to the four corners of the bounding box
    double dist[2][2] = { {INFINITY, INFINITY}, {INFINITY, INFINITY} };

    for (auto v_it = omMesh.vertices_begin();
              v_it != omMesh.vertices_end();  ++v_it)
    {
        const auto& pt = omMesh.point(v_it.handle());

        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            double d = std::hypot(pt[0] - box[0][i], pt[1] - box[1][j]);
            dist[i][j] = std::min(d, dist[i][j]);
        }
    }

    omMesh.request_vertex_status();

    // lock vertices closest to the corners
    // (note: multiple vertices can have the same minimum distance to a corner)
    for (auto v_it = omMesh.vertices_begin();
              v_it != omMesh.vertices_end();  ++v_it)
    {
        const auto& pt = omMesh.point(v_it.handle());

        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            double d = std::hypot(pt[0] - box[0][i], pt[1] - box[1][j]);
            if (std::abs(d - dist[i][j]) < 1e-12)
            {
                auto handle = v_it.handle();
                LOG(info1) << "Locking corner vertex " << handle.idx();
                omMesh.status(handle).set_locked(true);
            }
        }
    }
}

math::Extents2 prepareMesh(OMMesh &omMesh, const Mesh &mesh
                           , const SimplifyOptions &options)
{
    math::Extents2 me;

    // preprocess mesh - remove non-manifold edges
    if (options.check(SimplifyOption::RMNONMANIFOLDEDGES)) {
        auto pmesh(geometry::removeNonManifoldEdges(mesh));
        me = toOpenMesh(*pmesh, omMesh);
    } else{
        me = toOpenMesh(mesh, omMesh);
    }
    omMesh.update_normals();

    // lock the corner and/or border vertices based on flags
    if (options.check(SimplifyOption::CORNERS)) {
        lockCorners(omMesh);
    }
    if (options.check(SimplifyOption::INNERBORDER)) {
        lockBorder(omMesh, true, false);
    }
    if (options.check(SimplifyOption::OUTERBORDER)) {
        lockBorder(omMesh, false, true);
    }

    // return accumulated extents
    return me;
}

void prepareDecimator(Decimator &decimator
                      , const SimplifyOptions &options)
{
    if (options.alternativeVertices()) {
        if (options.concaveVertexModifier()) {
            LOGTHROW (err3, std::runtime_error) <<
                "Concave/Convex modification is not "
                "implemented for alternative verteces";
        }
        HModQuadricHybrid hModQuadricHybrid;
        decimator.add(hModQuadricHybrid);
        decimator.module(hModQuadricHybrid)
            .setAlternativeVertices(options.alternativeVertices());

        if (options.maxError()) {
            decimator.module(hModQuadricHybrid)
                .set_max_err(*options.maxError(), false);
        }
    } else if (options.concaveVertexModifier()){
        // collapse priority based on vertex error quadric
        // adjusted by convexity of the vertex
        HModQuadricConvex hModQuadricConvex;
        decimator.add(hModQuadricConvex);
        decimator.module(hModQuadricConvex).
            set_concave_vertex_modifier(*options.concaveVertexModifier());
        if (options.maxError()) {
            decimator.module(hModQuadricConvex)
                .set_max_err(*options.maxError(), false);
        }
    } else {
        // collapse priority based on vertex error quadric
        HModQuadric hModQuadric;
        decimator.add(hModQuadric);
        if (options.maxError()) {
            decimator.module(hModQuadric)
                .set_max_err(*options.maxError(), false);
        }
    }

    // apply normal flipping prevention (if configured)
    if (options.check(SimplifyOption::PREVENTFACEFLIP)) {
        HModNormalFlippingT hModNormalFlippingT;
        decimator.add(hModNormalFlippingT);
        decimator.module(hModNormalFlippingT).set_max_normal_deviation(90);
    }

    // apply aspect ratio (if configured)
    if (options.minAspectRatio()) {
        HModAspectRatioT h;
        decimator.add(h);
        decimator.module(h).set_aspect_ratio(*options.minAspectRatio());
    }

    // apply max edge lenght (if configured)
    if (options.maxEdgeLength()) {
        HModEdgeLengthT h;
        decimator.add(h);
        decimator.module(h).set_edge_length(*options.maxEdgeLength());
    }
}

} // namespace

Mesh::pointer simplify(const Mesh &mesh, int faceCount
                       , const SimplifyOptions &options)
{
    OMMesh omMesh;
    prepareMesh(omMesh, mesh, options);

    Decimator decimator(omMesh);
    prepareDecimator(decimator, options);
    decimator.initialize();

    decimator.decimate_to_faces(0, faceCount);
    omMesh.garbage_collection();

    auto newMesh(std::make_shared<geometry::Mesh>());
    fromOpenMesh(omMesh, *newMesh);
    return newMesh;
}

Mesh::pointer simplifyToError(const Mesh &mesh, double maxErr
                            , const SimplifyOptions &options)
{
    OMMesh omMesh;

    prepareMesh(omMesh, mesh, options);

    Decimator decimator(omMesh);
    HModQuadric hModQuadric; // collapse priority based on vertex error quadric
    decimator.add(hModQuadric);
    decimator.module( hModQuadric ).set_max_err(maxErr, false);
    decimator.initialize();

    decimator.decimate_to_faces(0, 0);
    omMesh.garbage_collection();

    auto newMesh(std::make_shared<geometry::Mesh>());
    fromOpenMesh(omMesh, *newMesh);
    return newMesh;
}

namespace {

inline std::size_t gridIndexImpl(const math::Size2_<long> &size
                                 , const math::Point2 origin
                                 , const math::Size2f tileSize
                                 , double x, double y)
{
    long xx((x - origin(0)) / tileSize.width);
    long yy((y - origin(1)) / tileSize.height);

    // clamp to grid (values can be outside of grid in case non-integer
    // tileSize
    xx = math::clamp(xx, long(0), long(size.width - 1));
    yy = math::clamp(yy, long(0), long(size.height -1));

    return xx + yy * size.width;
}

struct Tiling {
    math::Size2_<long> size;
    math::Point2 origin;
    math::Size2f tileSize;
    FacesPerCell::Functor facesPerCell;

    Tiling() : tileSize() {}

    Tiling(const math::Size2_<long> &size, const math::Point2 &origin
           , const math::Size2f &tileSize
           , const FacesPerCell::Functor &facesPerCell)
        : size(size), origin(origin), tileSize(tileSize)
        , facesPerCell(facesPerCell)
    {}

    inline std::size_t gridIndex(double x, double y) const {
        return gridIndexImpl(size, origin, tileSize, x, y);
    };

    inline math::Point2 tileLowerLeft(int index){
        int y = index / size.width;
        int x = index % size.width;

        return math::Point2(x * tileSize.width + origin(0)
                            , y * tileSize.height + origin(1));
    }

    inline std::vector<std::size_t> intersectingCells(math::Point2 triangle[3])
        const
    {
        std::vector<std::size_t> result;

        for(uint x = 0; x < size.width; ++x){
            for(uint y = 0; y < size.height; ++y){
                math::Point2 offset(x * tileSize.width
                                    , y * tileSize.height);
                math::Point2 ll = origin + offset;
                math::Point2 ur = ll + math::Point2
                    (tileSize.width, tileSize.height);

                if(math::triangleRectangleCollision(triangle,ll,ur)){
                    result.push_back(x + y * size.width);
                }
            }
        }
        return result;
    };

    /** Calculate maximum face count in all cells of grid. If real value is
     *  lower than maximum, set maximum to real value.
     */
    std::vector<std::size_t> getMax(const std::vector<std::size_t> &current)
        const
    {
        FacesPerCell fpc(size, origin, tileSize);
        facesPerCell(fpc);

        // merge minimum values from current and expected face counts
        auto out(fpc.get());
        auto icurrent(current.begin());
        for (auto &value : out) {
            value = std::min(value, *icurrent++);
        }

        // done
        return out;
    }

    static OpenMesh::MPropHandleT<Tiling> Property;
};

OpenMesh::MPropHandleT<Tiling> Tiling::Property;

template <typename MeshType>
class ModGrid : public OpenMesh::Decimater::ModBaseT<MeshType>
{
public:
    DECIMATING_MODULE(ModGrid, MeshType, Grid);

    ModGrid(MeshType &mesh)
        : Base(mesh, true) // true -> binary module!
        , mesh_(mesh), tiling_()
    {
        mesh_.add_property(cells_);
    }

    ~ModGrid() {
        mesh_.remove_property(cells_);
        //stat();
    }

    virtual void initialize() UTILITY_OVERRIDE {
        auto &prop(mesh_.mproperty(Tiling::Property));
        if (!prop.n_elements()) {
            LOGTHROW(err2, std::runtime_error)
                << "Tiling property not set at mesh.";
        }
        tiling_ = &prop[0];

        calculateGrid();
    }

    virtual float collapse_priority(const CollapseInfo &ci) UTILITY_OVERRIDE {
        // can be both faces removed?
        if (canRemove(ci.fl) && canRemove(ci.fr)) {
            return Base::LEGAL_COLLAPSE;
        }

        return Base::ILLEGAL_COLLAPSE;
    }

    virtual void preprocess_collapse(const CollapseInfo &ci) UTILITY_OVERRIDE {
        removed(ci.fl);
        removed(ci.fr);
    }

    virtual void postprocess_collapse(const CollapseInfo &ci) UTILITY_OVERRIDE
    {
        // traverse all remaining faces of non-removed vertex
        for (auto fi(mesh_.vf_begin(ci.v1)), fe(mesh_.vf_end(ci.v1));
             fi != fe; ++fi)
        {
            const auto f(fi.handle());
            update(f, cell(f), barycenterCell(f));
        }
    }

    std::size_t desiredCount() const {
        return std::accumulate(max_.begin(), max_.end(), std::size_t(0));
    }

private:
    typedef typename MeshType::FaceHandle FaceHandle;
    typedef typename MeshType::Scalar Scalar;

    void stat() const {
        auto ioriginal(original_.begin());
        auto icurrent(current_.begin());
        for (auto max : max_) {
            LOG(info2) << *ioriginal++ << " -> " << *icurrent++
                       << " (" << max << ")";
        }
    }

    void calculateGrid() {
        current_.clear();
        current_.resize(area(tiling_->size));

        for (const auto &f : mesh_.faces()) {
            auto cellIndex(barycenterCell(f));

            // assign cell index to face
            cell(f) = cellIndex;

            // one more face at given cell
            ++current_[cellIndex];
        }

        // get face count limits
        max_ = tiling_->getMax(current_);

        // for statistics
        original_.assign(current_.begin(), current_.end());
    }

    inline std::size_t barycenterCell(FaceHandle f) const {
        // calculate barycenter of face (i.e. average x and y coordinate)
        math::Point2 p;
        int valence(0);
        for (auto vi(mesh_.cfv_iter(f)); vi; ++vi) {
            const auto &v(mesh_.point(vi.handle()));
            p(0) += v[0];
            p(1) += v[1];
            ++valence;
        }

        p(0) /= valence;
        p(1) /= valence;

        return tiling_->gridIndex(p(0), p(1));
    }

    void removed(FaceHandle fh) {
        if (fh.is_valid()) {
            --current_[cell(fh)];
        }
    }

    void update(FaceHandle fh, std::size_t oldCell, std::size_t newCell) {
        if (!fh.is_valid() || (oldCell == newCell)) {
            return;
        }

        // new cell index, update
        --current_[oldCell];
        ++current_[newCell];
        cell(fh) = newCell;
    }

    inline bool canRemove(FaceHandle fh) const {
        if (!fh.is_valid()) { return true; }

        auto index(cell(fh));
        return (current_[index] > max_[index]);
    }

    inline Scalar& cell(FaceHandle fh) {
        return mesh_.property(cells_, fh);
    }

    inline const Scalar& cell(FaceHandle fh) const {
        return mesh_.property(cells_, fh);
    }

    MeshType &mesh_;
    Tiling *tiling_;
    std::vector<std::size_t> current_;
    std::vector<std::size_t> max_;
    std::vector<std::size_t> original_;
    OpenMesh::FPropHandleT<Scalar> cells_;
};


 /*
 * OpenMesh module for decimation mesh in grid.
 * Each cell has defined maximal face count and
 * once this face count is reached, incident faces
 * can't be decimated.
 */
template <typename MeshType>
class ModGrid2 : public OpenMesh::Decimater::ModBaseT<MeshType>
{
public:
    DECIMATING_MODULE(ModGrid2, MeshType, Grid);

    ModGrid2(MeshType &mesh)
        : Base(mesh, true) // true -> binary module!
        , mesh_(mesh), tiling_()
    {
        mesh_.add_property(cells_);
    }

    ~ModGrid2() {
        mesh_.remove_property(cells_);
        //stat();
    }

    virtual void initialize() UTILITY_OVERRIDE {
        auto &prop(mesh_.mproperty(Tiling::Property));
        if (!prop.n_elements()) {
            LOGTHROW(err2, std::runtime_error)
                << "Tiling property not set at mesh.";
        }
        tiling_ = &prop[0];

        calculateGrid();
    }

    virtual float collapse_priority(const CollapseInfo &ci) UTILITY_OVERRIDE {
        // can be both faces removed?
        if (canRemove(ci.fl) && canRemove(ci.fr)) {
            return Base::LEGAL_COLLAPSE;
        }

        return Base::ILLEGAL_COLLAPSE;
    }

    virtual void preprocess_collapse(const CollapseInfo &ci) UTILITY_OVERRIDE {
        removed(ci.fl);
        removed(ci.fr);
    }

    virtual void postprocess_collapse(const CollapseInfo &ci) UTILITY_OVERRIDE
    {
        // traverse all remaining faces of non-removed vertex
        for (auto fi(mesh_.vf_begin(ci.v1)), fe(mesh_.vf_end(ci.v1));
             fi != fe; ++fi)
        {
            const auto f(fi.handle());
            update(f);
        }
    }

    std::size_t desiredCount() const {
        return std::accumulate(max_.begin(), max_.end(), std::size_t(0));
    }

private:
    typedef typename MeshType::FaceHandle FaceHandle;
    typedef typename MeshType::Scalar Scalar;
    typedef typename std::vector<std::size_t> CellList;
    typedef typename std::vector<std::set<FaceHandle>> CellFaces;

    void setCellVerticesAsFeatures(int cellIndex) const{
        for(const auto &face : cellFaces_[cellIndex]){
            for (auto vi(mesh_.cfv_iter(face)); vi; ++vi) {
                mesh_.status(vi.handle()).set_feature(true);
            }
        }
    }

    void printCells() const {
        LOG(info2) << "Cells count:";
        for (uint i=0; i<current_.size(); ++i) {
            LOG(info2) << tiling_->tileLowerLeft(i)<<" - "<<current_[i];
        }
    }

    void stat() const {
        auto ioriginal(original_.begin());
        auto icurrent(current_.begin());
        for (auto max : max_) {
            LOG(info2) << *ioriginal++ << " -> " << *icurrent++
                       << " (" << max << ")";
        }
        printCells();
    }

    void calculateGrid() {
        current_.clear();
        current_.resize(area(tiling_->size));
        cellFaces_.clear();
        cellFaces_.resize(area(tiling_->size));

        for (const auto &f : mesh_.faces()) {
            cellList(f).clear();
            update(f);
        }

        // get face count limits
        max_ = tiling_->getMax(current_);

        // for statistics
        original_.assign(current_.begin(), current_.end());
    }

    void removed(FaceHandle fh) {
        if (fh.is_valid()) {
            for(auto &cell : cellList(fh)){
                --current_[cell];
                cellFaces_[cell].erase(fh);
            }
        }
    }

    void update(FaceHandle fh) {
        if (!fh.is_valid()) {
            return;
        }

        // decrement facecount in old cells
        for(auto &cell : cellList(fh)) {
            --current_[cell];
            cellFaces_[cell].erase(fh);
        }
        // recalculate cells
        calculateCells(fh);

        // increment facecount in new cells
        for(auto &cell : cellList(fh)) {
            ++current_[cell];
            cellFaces_[cell].insert(fh);
        }
    }

    inline void calculateCells(FaceHandle fh){
        math::Point2 vertices[3];
        uint vc = 0;
        for (auto vi(mesh_.cfv_iter(fh)); vi; ++vi) {
            const auto &v(mesh_.point(vi.handle()));
            vertices[vc](0) =  v[0];
            vertices[vc](1) =  v[1];
            vc++;
        }
        cellList(fh) = tiling_->intersectingCells(vertices);
    }

    inline bool canRemove(FaceHandle fh) const {
        if (!fh.is_valid()) { return true; }
        for(auto &cell : cellList(fh)){
            if(current_[cell] <= max_[cell]){
                setCellVerticesAsFeatures(cell);
                return false;
            }
        }
        return true;
    }

    inline CellList& cellList(FaceHandle fh) {
        return mesh_.property(cells_, fh);
    }

    inline const CellList& cellList(FaceHandle fh) const {
        return mesh_.property(cells_, fh);
    }

    MeshType &mesh_;
    Tiling *tiling_;
    std::vector<std::size_t> original_;
    std::vector<std::size_t> current_;
    std::vector<std::size_t> max_;
    OpenMesh::FPropHandleT<CellList> cells_;
    CellFaces cellFaces_;
};


} // namespace

/** Converts coordinate to grid defined by reference point and cell size.
 *  Coordinates not on grid are rounded down.
 */
double gridExtentsDown(double value, double origin, double size)
{
    return std::floor((value - origin) / size) * size + origin;
}

/** Converts coordinate to grid defined by reference point and cell size.
 *  Coordinates not on grid are rounded up.
 */
double gridExtentsUp(double value, double origin, double size)
{
    return std::ceil((value - origin) / size) * size + origin;
}

math::Extents2 gridExtents(const math::Extents2 &extents
                           , const math::Point2 &alignment
                           , const math::Size2f &cellSize)
{
    return {
        gridExtentsDown(extents.ll(0), alignment(0), cellSize.width)
        , gridExtentsDown(extents.ll(1), alignment(1), cellSize.height)
        , gridExtentsUp(extents.ur(0), alignment(0), cellSize.width)
        , gridExtentsUp(extents.ur(1), alignment(1), cellSize.height)
    };
}


std::size_t FacesPerCell::cellIndex(double x, double y) const
{
    return gridIndexImpl(gridSize_, origin_, cellSize_, x, y);
}

math::Extents2 FacesPerCell::cellExtents(std::size_t index) const
{
    // convert linear index to grid index
    auto i(index % gridSize_.width);
    auto j(index / gridSize_.width);

    // return extents
    return math::Extents2(origin_(0) + i * cellSize_.width
                          , origin_(1) + j * cellSize_.height
                          , origin_(0) + (i + 1) * cellSize_.width
                          , origin_(1) + (j + 1) * cellSize_.height);
}

Mesh::pointer simplifyInGrid(const Mesh &mesh, const math::Point2 &alignment
                             , const math::Size2f &cellSize
                             , const FacesPerCell::Functor &facesPerCell
                             , const SimplifyOptions &options)
{
    LOG(info2) << "[simplify] alignment: "
               << std::setprecision(15) << alignment;
    LOG(info2) << "[simplify] cellSize: " << cellSize;

    // convert to openmesh structure
    OMMesh omMesh;
    math::Extents2 me(prepareMesh(omMesh, mesh, options));

    // calculate grid extents
    const auto ge(gridExtents(me, alignment, cellSize));

    // grid origin
    const auto gorigin(ge.ll);

    // grid size
    const math::Size2_<long> gsize
        (long(std::round(((ge.ur(0) - ge.ll(0)) / cellSize.width)))
         , long(std::round((ge.ur(1) - ge.ll(1)) / cellSize.width)));

    LOG(info2) << "[simplify] mesh extents: " << me;
    LOG(info2) << "[simplify] gridded mesh extents: "
               << std::setprecision(15) << ge;
    LOG(info2) << "[simplify] grid size: " << gsize;
    LOG(info2) << "[simplify] grid origin: " << gorigin;

    // create and add tiling as a mesh property to mesh
    omMesh.add_property(Tiling::Property);
    omMesh.mproperty(Tiling::Property).push_back();
    omMesh.mproperty(Tiling::Property)[0]
        = Tiling(gsize, gorigin, cellSize, facesPerCell);

    // create and prepare decimator
    Decimator decimator(omMesh);
    prepareDecimator(decimator, options);

    // add decimator grid module
    ModGrid<OMMesh>::Handle hModGrid;
    decimator.add(hModGrid);

    // initialized
    decimator.initialize();

    auto fc(decimator.module(hModGrid).desiredCount());
    LOG(info2) << "[simplify] Simplifying mesh to " << fc << " faces.";
    // simplifying to 0 faces, simplification will stop based on grid constraints
    decimator.decimate_to_faces(0, 0);

    omMesh.garbage_collection();

    auto newMesh(std::make_shared<geometry::Mesh>());
    fromOpenMesh(omMesh, *newMesh);
    return newMesh;
}


} // namespace geometry

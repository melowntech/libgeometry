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

#include "./meshop.hpp"

namespace geometry {

/* code moved from the window-mesh utility */

namespace {

typedef OpenMesh::TriMesh_ArrayKernelT<> OMMesh;
typedef OpenMesh::Decimater::DecimaterT<OMMesh> Decimator;
typedef OpenMesh::Decimater::ModQuadricT<OMMesh>::Handle HModQuadric;

/** Convert mesh to OpenMesh data structure and return mesh extents (2D bounding
 *  box)
 */
math::Extents2 toOpenMesh(const geometry::Mesh &imesh, OMMesh& omMesh)
{
    math::Extents2 e(math::InvalidExtents{});

    //preprocess mesh - remove non-manifold edges
    auto pmesh(removeNonManifoldEdges(imesh));
    auto & mesh(*pmesh);
    uint removedFaces = imesh.faces.size() - mesh.faces.size();
    if(removedFaces != 0){
        LOG( info2 )<< "Removed " <<removedFaces
                    << " faces (incidental to non-manifold edge)";
    }

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

Mesh::pointer simplify(const Mesh &mesh, int faceCount)
{
    OMMesh omMesh;
    toOpenMesh(mesh, omMesh);

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

Mesh::pointer simplifyToError(const Mesh &mesh, double maxErr)
{
    OMMesh omMesh;
    toOpenMesh(mesh, omMesh);

    // lock the corner vertices of the window to prevent simplifying the corners
    lockCorners(omMesh);

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

struct Tiling {
    math::Size2_<long> size;
    math::Point2 origin;
    double tileSize;
    FacesPerCell facesPerCell;

    Tiling() : tileSize() {}

    Tiling(const math::Size2_<long> &size, const math::Point2 &origin
           , double tileSize, const FacesPerCell &facesPerCell)
        : size(size), origin(origin), tileSize(tileSize)
        , facesPerCell(facesPerCell)
    {}

    inline std::size_t gridIndex(double x, double y) const {
        long xx((x - origin(0)) / tileSize);
        long yy((y - origin(1)) / tileSize);

#ifndef NDEBUG
        // sanity check (only in debug mode)
        if ((xx < 0) || (xx >= size.width)
            || (yy < 0) || (yy >= size.height))
        {
            LOGTHROW(err2, std::runtime_error)
                << "Grid index out of bounds!";
        }
#endif

        return xx + yy * size.width;
    };

    inline math::Point2 tileLowerLeft(int index){
        int y = index/size.width;
        int x = index%size.width;

        return math::Point2( x*tileSize+origin(0)
                           , y*tileSize+origin(1));
    }

    inline std::vector<std::size_t> intersectingCells(math::Point2 triangle[3]) const{
        std::vector<std::size_t> result;

        for(uint x=0; x<size.width; ++x){
            for(uint y=0; y<size.height; ++y){
                math::Point2 offset(x*tileSize, y*tileSize);
                math::Point2 ll = origin + offset;
                math::Point2 ur = ll+math::Point2(tileSize, tileSize);

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
    std::vector<std::size_t> getMax(const std::vector<std::size_t> &real)
        const
    {
        std::vector<std::size_t> max;

        auto ireal(real.begin());
        for (long j(0); j < size.height; ++j) {
            for (long i(0); i < size.width; ++i, ++ireal) {
                max.push_back
                    (std::min
                     (facesPerCell
                      (math::Extents2(origin(0) + i * tileSize
                                      , origin(1) + j * tileSize
                                      , origin(0) + (i + 1) * tileSize
                                      , origin(1) + (j + 1) * tileSize))
                      , *ireal));
            }
        }

        return max;
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

        //decrement facecount in old cells
        for(auto &cell : cellList(fh)){
            --current_[cell];
            cellFaces_[cell].erase(fh);
        }
        //recalculate cells
        calculateCells(fh);

        //increment facecount in new cells
        for(auto &cell : cellList(fh)){
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


/** Converts coordinate to grid defined by reference point and cell size.
 *  Coordinates not on grid are rounded up.
 */
double gridExtentsUp(double value, double origin, double size)
{
    return std::ceil((value - origin) / size) * size + origin;
}

/** Converts coordinate to grid defined by reference point and cell size.
 *  Coordinates not on grid are rounded down.
 */
double gridExtentsDown(double value, double origin, double size)
{
    return std::floor((value - origin) / size) * size + origin;
}

math::Extents2 gridExtents(const math::Extents2 &extents
                           , const math::Point2 &alignment
                           , double size)
{
    return {
        gridExtentsDown(extents.ll(0), alignment(0), size)
        , gridExtentsDown(extents.ll(1), alignment(1), size)
        , gridExtentsUp(extents.ur(0), alignment(0), size)
        , gridExtentsUp(extents.ur(1), alignment(1), size)
    };
}

} // namespace

Mesh::pointer simplifyInGrid(const Mesh &mesh, const math::Point2 &alignment
                             , double cellSize
                             , const FacesPerCell &facesPerCell)
{
    LOG(info2) << "[simplify] alignment: "
               << std::setprecision(15) << alignment;
    LOG(info2) << "[simplify] cellSize: " << cellSize;

    // convert tom openmehs structure
    OMMesh omMesh;
    const auto me(toOpenMesh(mesh, omMesh));

    // calculate grid extents
    const auto ge(gridExtents(me, alignment, cellSize));

    // grid origin
    const auto gorigin(ge.ll);

    // grid size
    const math::Size2_<long> gsize
        (long((ge.ur(0) - ge.ll(0)) / cellSize)
         , long((ge.ur(1) - ge.ll(1)) / cellSize));

    LOG(info2) << "[simplify] mesh extents: " << me;
    LOG(info2) << "[simplify] gridded mesh extents: "
               << std::setprecision(15) << ge;
    LOG(info2) << "[simplify] grid size: " << gsize;
    LOG(info2) << "[simplify] grid origin: " << gorigin;

    // lock the corner vertices of the window to prevent simplifying the corners
    lockCorners(omMesh);

    // create and add tiling as a mesh property to mesh
    omMesh.add_property(Tiling::Property);
    omMesh.mproperty(Tiling::Property).push_back();
    omMesh.mproperty(Tiling::Property)[0]
        = Tiling(gsize, gorigin, cellSize, facesPerCell);

    Decimator decimator(omMesh);

    // collapse priority based on vertex error quadric
    HModQuadric hModQuadric;
    decimator.add(hModQuadric);

    ModGrid2<OMMesh>::Handle hModGrid;
    decimator.add(hModGrid);

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

/**
 *  @file nnquadtree.hpp
 *  @author Jakub Cerveny <jakub.cerveny@melown.com>
 *
 *  Quadtree for nearest neighbor search.
 */

#ifndef geometry_nnquadtree_hpp_included_
#define geometry_nnquadtree_hpp_included_

#include <limits>
#include <cmath>

namespace geometry {

/** NNQuadTree -- quadtree that stores 2D points and provides nearest
 *  neighbor search. Points are first inserted into the tree (one by one),
 *  then nearest neighbor queries can be issued.
 *
 *  PointType -- must define operator() for coordinate access (e.g.,
 *               math::Point2f)
 *  IdType    -- type for point IDs (e.g, int)
 *  NodeSize  -- maximum number of points stored in a leaf node
 */
template<typename PointType, typename IdType, int NodeSize = 4>
class NNQuadTree
{
public:

    /** Initialize the tree, min and max are the bounding rectangle where all
     *  points must lie. This is the rectangle that will be subdivided by the
     *  tree.
     */
    NNQuadTree(const PointType &min, const PointType &max)
        : root_(new Node), min_(min), max_(max)
    {}

    /** Insert a point into the tree. For future reference, an ID of the point
     *  needs to also be passed.
     */
    void insert(const PointType &point, IdType id)
    {
        root_->insert(point, id, min_, max_);
    }

    /** Find the nearest point to 'query' and return its ID. If the tree is
     *  empty, an invalid ID is returned (can be changed with 'invalidId').
     */
    IdType nearest(const PointType &query,
                   IdType invalidId = std::numeric_limits<IdType>::max())
    {
        IdType minId = invalidId;
        double minDist = std::numeric_limits<double>::max();
        root_->nearest(query, min_, max_, minDist, minId);
        return minId;
    }

    /** Find the nearest point to 'query' and return the distance between
     *  the two. If the tree is empty, maximum 'double' is returned.
     */
    double nearestDist(const PointType &query)
    {
        IdType minId = IdType();
        double minDist = std::numeric_limits<double>::max();
        root_->nearest(query, min_, max_, minDist, minId);
        return minDist;
    }

    ~NNQuadTree() { delete root_; }

protected:

    struct Node
    {
        union {
            Node* ch[4];
            struct {
                PointType pt[NodeSize];
                IdType id[NodeSize];
            } data;
        };
        int size; // >= 0 if leaf, < 0 if internal node

        Node() : size(0) {}

        void insert(const PointType &point, IdType id,
                    const PointType &min, const PointType &max);

        void nearest(const PointType &query,
                     const PointType &min, const PointType &max,
                     double &minDist, IdType &minId);

        ~Node() {
            if (size < 0) {
                for (int i = 0; i < 4; i++) {
                    delete ch[i];
                }
            }
        }
    };

    Node* root_;
    PointType min_, max_; // root extents

    // TODO: a fast allocator for Nodes
    // TODO: user shouldn't need to specify initial bounding box
    // TODO: store bounding box in slack space of branch nodes?
};


// implementation

template<typename PointType, typename IdType, int NodeSize>
void NNQuadTree<PointType, IdType, NodeSize>::Node::insert(
        const PointType &point, IdType id,
        const PointType &min, const PointType &max)
{
    if (size < 0) // branch node
    {
        double x0 = min(0), x2 = max(0), x1 = (x0 + x2)*0.5;
        double y0 = min(1), y2 = max(1), y1 = (y0 + y2)*0.5;

        //    y2 *-------*-------* max
        //       |       |       |
        //       |   2   |   3   |
        //       |       |       |
        //    y1 *-------*-------*
        //       |       |       |
        //       |   0   |   1   |
        //       |       |       |
        //    y0 *-------*-------*
        //       x0      x1      x2

        if (point(0) < x1) {
            if (point(1) < y1) {
                ch[0]->insert(point, id, min, PointType(x1, y1));
            }
            else {
                ch[2]->insert(point, id, PointType(x0, y1), PointType(x1, y2));
            }
        }
        else {
            if (point(1) < y1) {
                ch[1]->insert(point, id, PointType(x1, y0), PointType(x2, y1));
            }
            else {
                ch[3]->insert(point, id, PointType(x1, y1), max);
            }
        }
    }
    else if (size < NodeSize) // leaf, not full
    {
        data.pt[size] = point;
        data.id[size] = id;
        size++;
    }
    else // full leaf, we need to split it
    {
        PointType tmp_pt[NodeSize];
        IdType tmp_id[NodeSize];

        // back up node content
        for (int i = 0; i < NodeSize; i++) {
            tmp_pt[i] = data.pt[i];
            tmp_id[i] = data.id[i];
        }

        // turn this leaf into an internal node, create new empty children
        size = -1;
        for (int i = 0; i < 4; i++) {
            ch[i] = new Node;
        }

        // reinsert old points
        for (int i = 0; i < NodeSize; i++) {
            insert(tmp_pt[i], tmp_id[i], min, max);
        }

        // finally, insert the new one
        insert(point, id, min, max);
    }
}


template<typename PointType, typename IdType, int NodeSize>
void NNQuadTree<PointType, IdType, NodeSize>::Node::nearest(
        const PointType &query,
        const PointType &min, const PointType &max,
        double &minDist, IdType &minId)
{
    // quickly reject subtrees that cannot contain a closer point
    if (min(0) - query(0) > minDist ||
        min(1) - query(1) > minDist ||
        query(0) - max(0) > minDist ||
        query(1) - max(1) > minDist)
    {
        return;
    }

    if (size < 0) // branch node, recurse to subtrees (closest first)
    {
        double x[3] = { min(0), (min(0) + max(0))*0.5, max(0) };
        double y[3] = { min(1), (min(1) + max(1))*0.5, max(1) };

        int qx = (query(0) < x[1]) ? 0 : 1;
        int qy = (query(1) < y[1]) ? 0 : 1;

        int order[4] = { 2*qy + qx,
                         2*qy + (1-qx),
                         2*(1-qy) + qx,
                         2*(1-qy) + (1-qx) };

        for (int i = 0; i < 4; i++)
        {
            int q = order[i];
            ch[q]->nearest(query,
                           PointType(x[q & 1], y[q >> 1]),
                           PointType(x[(q & 1) + 1], y[(q >> 1) + 1]),
                           minDist, minId);
        }
    }
    else // leaf node
    {
        for (int i = 0; i < size; i++)
        {
            double d2 = math::sqr(data.pt[i](0) - query(0)) +
                        math::sqr(data.pt[i](1) - query(1));

            if (d2 < math::sqr(minDist)) {
                minDist = std::sqrt(d2);
                minId = data.id[i];
            }
        }
    }
}

} // namespace geometry

#endif // geometry_nnquadtree_hpp_included_

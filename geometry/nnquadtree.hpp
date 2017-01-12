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

/** NNQuadTree -- quadtree that stores 2D points and provides a nearest
 *  neighbor search. Points are first inserted into the tree (one by one),
 *  then nearest neighbor queries can be issued.
 *
 *  PointType must define operator() for coordinate access (e.g., math::Point2f)
 *  IdType   -- type for point IDs (e.g, int)
 *  NodeSize -- maximum number of points stored in a leaf node
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
                ch[0]->insert(point, id, min, {x1, y1});
            }
            else {
                ch[2]->insert(point, id, {x0, y1}, {x1, y2});
            }
        }
        else {
            if (point(1) < y1) {
                ch[1]->insert(point, id, {x1, y0}, {x2, y1});
            }
            else {
                ch[3]->insert(point, id, {x1, y1}, max);
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
        // back up node content
        PointType tmp_pt[NodeSize];
        IdType tmp_id[NodeSize];
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
    if (size < 0) // branch node
    {
        double x0 = min(0), x2 = max(0), x1 = (x0 + x2)*0.5;
        double y0 = min(1), y2 = max(1), y1 = (y0 + y2)*0.5;

        if (query(0) < x1) {
            if (query(1) < y1) {
                ch[0]->nearest(query, min, {x1, y1}, minDist, minId);

                double xgap = math::sqr(x1 - query(0));
                double ygap = math::sqr(y1 - query(1));

                if (xgap < minDist) {
                    ch[1]->nearest(query, {x1, y0}, {x2, y1}, minDist, minId);
                }
                if (ygap < minDist) {
                    ch[2]->nearest(query, {x0, y1}, {x1, y2}, minDist, minId);
                }
                if (xgap < minDist && ygap < minDist) {
                    ch[3]->nearest(query, {x1, y1}, max, minDist, minId);
                }
            }
            else {
                ch[2]->nearest(query, {x0, y1}, {x1, y2}, minDist, minId);

                double xgap = math::sqr(x1 - query(0));
                double ygap = math::sqr(query(1) - y1);

                if (xgap < minDist) {
                    ch[3]->nearest(query, {x1, y1}, max, minDist, minId);
                }
                if (ygap < minDist) {
                    ch[0]->nearest(query, min, {x1, y1}, minDist, minId);
                }
                if (xgap < minDist && ygap < minDist) {
                    ch[1]->nearest(query, {x1, y0}, {x2, y1}, minDist, minId);
                }
            }
        }
        else {
            if (query(1) < y1) {
                ch[1]->nearest(query, {x1, y0}, {x2, y1}, minDist, minId);

                double xgap = math::sqr(query(0) - x1);
                double ygap = math::sqr(y1 - query(1));

                if (xgap < minDist) {
                    ch[0]->nearest(query, min, {x1, y1}, minDist, minId);
                }
                if (ygap < minDist) {
                    ch[3]->nearest(query, {x1, y1}, max, minDist, minId);
                }
                if (xgap < minDist && ygap < minDist) {
                    ch[2]->nearest(query, {x0, y1}, {x1, y2}, minDist, minId);
                }
            }
            else {
                ch[3]->nearest(query, {x1, y1}, max, minDist, minId);

                double xgap = math::sqr(query(0) - x1);
                double ygap = math::sqr(query(1) - y1);

                if (xgap < minDist) {
                    ch[2]->nearest(query, {x0, y1}, {x1, y2}, minDist, minId);
                }
                if (ygap < minDist) {
                    ch[1]->nearest(query, {x1, y0}, {x2, y1}, minDist, minId);
                }
                if (xgap < minDist && ygap < minDist) {
                    ch[0]->nearest(query, min, {x1, y1}, minDist, minId);
                }
            }
        }
    }
    else // leaf node
    {
        for (int i = 0; i < size; i++)
        {
            double d = math::sqr(data.pt[i](0) - query(0)) +
                       math::sqr(data.pt[i](1) - query(1));

            if (d < minDist) {
                minDist = d;
                minId = data.id[i];
            }
        }
    }
}

} // namespace geometry

#endif // geometry_nnquadtree_hpp_included_

/**
 *  @file kdtree.hpp
 *  @author Jakub Cerveny <jakub.cerveny@ext.citationtech.net>
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  The k-d tree data structure for fast nearest neighbor search.
 *
 *  2011-09-28 (vasek)    refactored
 */

#ifndef KDTREE_HPP_INCLUDED
#define KDTREE_HPP_INCLUDED

#include <vector>
#include <algorithm>

#include <boost/noncopyable.hpp>

namespace geometry {


/**
 *  KdTree -- implementation of a k-d tree. The k-d tree is a generalization of
 *  the binary search tree which enables fast nearest neighbor searches in
 *  more dimensions (k > 1). See http://en.wikipedia.org/wiki/K-d_tree for
 *  more details.
 *
 *  The type T is a point in space and K is the number of dimensions.
 *  The coordinates of T must be accessible with operator().
 *
 *  NB: Associated vector must outlive this tree!
 *  Associated vector is not modified by any means.
 */

namespace detail {

template<typename T, unsigned int K = 3>
class KdTreeNode : boost::noncopyable
{
public:
    typedef typename std::vector<T>::const_iterator iterator;

    typedef std::vector<iterator> Indirect;

    KdTreeNode(typename Indirect::iterator beg, typename Indirect::iterator end
               , int depth = 0)
        : sons({ nullptr, nullptr })
    {
        // trivial case - just one point
        if (beg + 1 >= end) {
            point = *beg;
            return;
        }

        // sort points by one of the coordinates
        unsigned int axis = depth % K;
        std::sort(beg, end, [axis](const iterator& a, const iterator& b)
                  { return (*a)(axis) < (*b)(axis); } );

        // the median will be the content (and boundary) of this node
        auto count(std::distance(beg, end));
        auto median(beg + count / 2);
        point = *median;

        // create two subtrees (points smaller and larger than the median)
        if (median > beg) {
            sons[0] = new KdTreeNode(beg, median, depth + 1);
        }
        if (median + 1 < end) {
            sons[1] = new KdTreeNode(median + 1, end, depth + 1);
        }
    }

    /**
     *  Finds the nearest neighbor of the point "query" (this is a general point
     *  which doesn't need to be in the tree). Returns the closest point from
     *  the tree and "dist2" is filled with the square of the distance.
     *  "axis" must be zero.
     */
    template<bool IgnoreEqual = false>
    iterator nearest(const T& query, double& dist2, unsigned int axis = 0)
        const
    {
        // calculate the squared distance between query and us
        T diff(query - *point);
        dist2 = inner_prod(diff, diff);

        // if a point already in the tree is queried we may want to ignore it
        if (IgnoreEqual && dist2 == 0.0) {
            dist2 = INFINITY;
        }

        // if this is a leaf, we're done
        iterator result = point;
        if (!sons[0] && !sons[1]) {
            return result;
        }

        // perpendicular distance to node boundary
        double perp = query(axis) - (*point)(axis);

        // change axis
        if (++axis >= K) {
            axis = 0;
        }

        double sub2;

        // search the near side of the tree first
        unsigned int side = (perp < 0) ? 0 : 1;
        if (sons[side]) {
            iterator it = sons[side]->nearest<IgnoreEqual>(query, sub2, axis);
            if (sub2 < dist2) {
                result = it;
                dist2 = sub2;

                // the already found distance may be so small that there's no
                // need to search the other sub-tree (we can't cross the boundary)
                if (sub2 < perp*perp)
                    return result;
            }
        }

        // we have to search the other side too
        side ^= 1;
        if (sons[side]) {
            iterator it = sons[side]->nearest<IgnoreEqual>(query, sub2, axis);
            if (sub2 < dist2) {
                result = it;
                dist2 = sub2;
            }
        }

        return result;
    }

    /**
     *  Returns all points that are within "radius" from "query".
     */
    template<bool IgnoreEqual = false>
    void range(const T& query, double radius,
               std::vector<T>& result, unsigned int axis = 0) const
    {
        // return point if it is within radius
        T diff(query - *point);
        double dist2 = inner_prod(diff, diff);
        if (dist2 <= radius*radius) {
            if (!IgnoreEqual || dist2 > 0.0)
                result.push_back(*point);
        }

        // perpendicular distance to node boundary
        double perp = query(axis) - (*point)(axis);

        // change axis
        if (++axis >= K) {
            axis = 0;
        }

        // recurse to sub-trees if they are within radius
        if (sons[0] && perp <= +radius) {
            sons[0]->range<IgnoreEqual>(query, radius, result, axis);
        }
        if (sons[1] && perp >= -radius) {
            sons[1]->range<IgnoreEqual>(query, radius, result, axis);
        }
    }

    ~KdTreeNode() {
        delete sons[0];
        delete sons[1];
    }

protected:
    static double inner_prod(const T &op1, const T &op2) {
        double retval( 0.0 );
        for (unsigned int i = 0; i < K; i++) {
            retval += op1(i) * op2(i);
        }
        return retval;
    };


    iterator point;
    KdTreeNode *sons[2];

    /* TODO: since the tree is perfectly balanced its structure is known
       implicitly and we could operate just on the point array with no
       additional nodes and pointers */
};

} // namespace detail

template<typename T, unsigned int K = 3>
class KdTree : boost::noncopyable
{
public:
    typedef typename std::vector<T>::const_iterator iterator;

    KdTree(iterator beg, iterator end)
        : root_(nullptr)
    {
        typename detail::KdTreeNode<T, K>::Indirect indirect;
        indirect.reserve(std::distance(beg, end));
        for ( ; beg != end; ++beg) {
            indirect.push_back(beg);
        }

        root_ = new detail::KdTreeNode<T, K>(indirect.begin(), indirect.end());
    }

    template<bool IgnoreEqual = false>
    iterator nearest(const T& query, double& dist2) const
    {
        return root_->nearest(query, dist2);
    }

    template<bool IgnoreEqual = false>
    std::vector<T> range(const T& query, double radius) const {
        std::vector<T> result;
        root_->range<IgnoreEqual>(query, radius, result);
        return result;
    }

    template<bool IgnoreEqual = false>
    std::pair<iterator, double> nearest(const T& query) const
    {
        std::pair<iterator, double> res;
        res.first = root_->nearest<IgnoreEqual>(query, res.second);
        return res;
    }

    template<bool IgnoreEqual = false>
    void range(const T& query, double radius, std::vector<T>& result) const {
        return root_->range<IgnoreEqual>(query, radius, result);
    }

    ~KdTree() {
        delete root_;
    }

private:
    detail::KdTreeNode<T, K> *root_;
};

} // namespace geometry

#endif

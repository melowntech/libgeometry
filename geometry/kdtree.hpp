/**
 *  @file kdtree.hpp
 *  @author Jakub Cerveny <jakub.cerveny@ext.citationtech.net>
 *
 *  The k-d tree data structure for fast nearest neighbor search.
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
 */

template<typename T, int K = 3>
class KdTree : boost::noncopyable
{
public:

    typedef typename std::vector<T>::iterator iterator;

    /**
     *  Builds the tree from a vector of points. You need to pass the beginning
     *  and end iterators of the vector. The vector is not read-only -- it will
     *  be modified (sorted) by this constructor. "depth" must be zero.
     */
    KdTree(iterator beg, iterator end, int depth = 0)
        : sons({ nullptr, nullptr })
    {
        // trivial case - just one point
        if (beg+1 >= end)
        {
            point = beg;
            return;
        }

        // sort points by one of the coordinates
        int axis = depth % K;
        std::sort(beg, end, [axis](const T& a, const T& b)
                                  { return a(axis) < b(axis); } );

        // the median will be the content (and boundary) of this node
        int count = end - beg;
        iterator median = beg + count/2;
        point = median;

        // create two subtrees (points smaller and larger than the median)
        if (median > beg) sons[0] = new KdTree(beg, median, depth+1);
        if (median+1 < end) sons[1] = new KdTree(median+1, end, depth+1);
    }


    /**
     *  Finds the nearest neighbor of the point "query" (this is a general point
     *  which doesn't need to be in the tree). Returns the closest point from
     *  the tree and "dist2" is filled with the square of the distance.
     *  "axis" must be zero.
     */
    template<int IgnoreEqual = false>
    iterator nearest(const T& query, double& dist2, int axis = 0) const
    {
        // calculate the squared distance between query and us
        T diff(query - *point);
        dist2 = inner_prod(diff, diff);

        // if a point already in the tree is queried we may want to ignore it
        if (IgnoreEqual && dist2 == 0.0)
            dist2 = INFINITY;

        // if this is a leaf, we're done
        iterator result = point;
        if (!sons[0] && !sons[1])
            return result;

        // perpendicular distance to node boundary
        double perp = query(axis) - (*point)(axis);

        // change axis
        if (++axis >= K) axis = 0;

        double sub2;

        // search the near side of the tree first
        int side = (perp < 0) ? 0 : 1;
        if (sons[side])
        {
            iterator it = sons[side]->nearest<IgnoreEqual>(query, sub2, axis);
            if (sub2 < dist2)
            {
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
        if (sons[side])
        {
            iterator it = sons[side]->nearest<IgnoreEqual>(query, sub2, axis);
            if (sub2 < dist2)
            {
                result = it;
                dist2 = sub2;
            }
        }

        return result;
    }

    template<int IgnoreEqual = false>
    std::pair<iterator, double> nearest(const T& query) const
    {
        std::pair<T, double> res;
        res.first = nearest<IgnoreEqual>(query, res.second);
        return res;
    }

    /**
     *  Returns all points that are within "radius" from "query".
     */
    template<int IgnoreEqual = false>
    void range(const T& query, double radius,
               std::vector<T>& result, int axis = 0) const
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
        if (++axis >= K) axis = 0;

        // recurse to sub-trees if they are within radius
        if (sons[0] && perp <= +radius) {
            sons[0]->range<IgnoreEqual>(query, radius, result, axis);
        }
        if (sons[1] && perp >= -radius) {
            sons[1]->range<IgnoreEqual>(query, radius, result, axis);
        }
    }

    ~KdTree()
    {
        delete sons[0];
        delete sons[1];
    }

protected:

    static double inner_prod( const T & op1, const T & op2 ) {
        
        double retval( 0.0 );
        for ( int i = 0; i < K; i++ )
            retval += op1(i) * op2(i);
        return retval;
    };
    

    iterator point;
    KdTree* sons[2];

    /* TODO: since the tree is perfectly balanced its structure is known
       implicitly and we could operate just on the point array with no
       additional nodes and pointers */
};


} // namespace geometry

#endif


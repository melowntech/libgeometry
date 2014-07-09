/**
 *  @file kdtree.hpp
 *  @author Jakub Cerveny <jakub.cerveny@ext.citationtech.net>
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  The k-d tree data structure for fast nearest neighbor search.
 *
 *  2012-09-28 (vasek)    refactored
 */

#ifndef KDTREE_HPP_INCLUDED
#define KDTREE_HPP_INCLUDED

#include <vector>
#include <algorithm>
#include <utility>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range/iterator_range.hpp>


namespace geometry {

/** Default coordination accessor.
 */
template <typename T> struct GetCoordinate;

/**
 *  KdTree -- implementation of a k-d tree. The k-d tree is a generalization of
 *  the binary search tree which enables fast nearest neighbor searches in
 *  more dimensions (k > 1). See http://en.wikipedia.org/wiki/K-d_tree for
 *  more details.
 *
 *  typename T:     point in space
 *  unsinged int K: number of dimensions
 *  typename G:     dimension value getter (defaults to T::operator())
 *  typaneme C:     associated container type (defaults to std::vector)
 *
 *  NB: Associated container must outlive this tree (tree holds iterators to it)
 *      Associated container is not modified by any means.
 */
template<typename T, unsigned int K = 3, typename G = GetCoordinate<T>
         , typename C = std::vector<T>>
class KdTree;

template<typename T, unsigned int K, typename G, typename C>
double samplingDelta(const KdTree<T, K, G, C> &tree
                     , double bulkThreshold = .5);

struct IntrusiveKdTree {};

namespace detail {

/**
 *  KdTreeNode: Actual tree node implementation.
 *
 *  typename T:     point in space
 *  unsinged int K: number of dimensions
 *  typename G:     dimension value getter (defaults to T::operator())
 *  typaneme C:     associated container type (defaults to std::vector)
 *
 *  NB: Associated iterators must not be invalidated/modified.
 *      Associated iterators are not modified by any means.
 */
template<typename T, unsigned int K, typename G, typename C>
class KdTreeNode : boost::noncopyable, G
{
public:
    typedef C container_type;
    typedef typename container_type::iterator mutable_iterator;
    typedef typename container_type::const_iterator iterator;

    typedef std::vector<iterator> Indirect;

    KdTreeNode(typename Indirect::iterator beg, typename Indirect::iterator end
               , int depth = 0)
        : sons { nullptr, nullptr }
    {
        // trivial case - just one point
        if (beg + 1 >= end) {
            point = *beg;
            return;
        }

        // sort points by one of the coordinates
        unsigned int axis = depth % K;
        std::sort(beg, end, [axis, this](const iterator& a, const iterator& b)
                  { return (this->G::get(*a, axis)
                            < this->G::get(*b, axis)); } );

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

    KdTreeNode(mutable_iterator beg, mutable_iterator end
               , const IntrusiveKdTree &intrusive, int depth = 0)
        : sons { nullptr, nullptr }
    {
        // trivial case - just one point
        if (beg + 1 >= end) {
            point = beg;
            return;
        }

        // sort points by one of the coordinates
        unsigned int axis = depth % K;
        std::sort(beg, end, [axis, this](const T &a, const T &b)
                  { return (this->G::get(a, axis)
                            < this->G::get(b, axis)); } );

        // the median will be the content (and boundary) of this node
        auto count(std::distance(beg, end));
        auto median(beg + count / 2);
        point = median;

        // create two subtrees (points smaller and larger than the median)
        if (median > beg) {
            sons[0] = new KdTreeNode(beg, median, intrusive, depth + 1);
        }
        if (median + 1 < end) {
            sons[1] = new KdTreeNode(median + 1, end, intrusive, depth + 1);
        }
    }

    /**
     *  Finds the nearest neighbor of the point "query" (this is a general point
     *  which doesn't need to be in the tree). Returns the closest point from
     *  the tree and "dist2" is filled with the square of the distance.
     *  "axis" must be zero.
     */
    template<bool IgnoreEqual>
    iterator nearest(const T& query, double& dist2, unsigned int axis = 0)
        const
    {
        // calculate the squared distance between query and us
        T diff(G::diff(query, *point));
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
        double perp = G::get(query, axis) - G::get(*point, axis);

        // change axis
        if (++axis >= K) {
            axis = 0;
        }

        double sub2;

        // search the near side of the tree first
        unsigned int side = (perp < 0) ? 0 : 1;
        if (sons[side]) {
            iterator it = sons[side]->template nearest<IgnoreEqual>(query, sub2, axis);
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
            iterator it = sons[side]->template nearest<IgnoreEqual>(query, sub2, axis);
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
    template<bool IgnoreEqual>
    void range(const T& query, double radius,
               std::vector<T>& result, unsigned int axis = 0) const
    {
        // return point if it is within radius
        T diff(G::diff(query, *point));
        double dist2 = inner_prod(diff, diff);
        if (dist2 <= radius*radius) {
            if (!IgnoreEqual || dist2 > 0.0)
                result.push_back(*point);
        }

        // perpendicular distance to node boundary
        double perp = G::get(query, axis) - G::get(*point, axis);

        // change axis
        if (++axis >= K) {
            axis = 0;
        }

        // recurse to sub-trees if they are within radius
        if (sons[0] && perp <= +radius) {
            sons[0]->template range<IgnoreEqual>(query, radius, result, axis);
        }
        if (sons[1] && perp >= -radius) {
            sons[1]->template range<IgnoreEqual>(query, radius, result, axis);
        }
    }

    /**
     *  Returns all points that are within "radius" from "query".
     *  The points are returned indirectly through iterators.
     */
    template<bool IgnoreEqual>
    void range(const T& query, double radius,
              std::vector<iterator>& result, unsigned int axis = 0) const
    {
       // return point if it is within radius
       T diff(G::diff(query, *point));
       double dist2 = inner_prod(diff, diff);
       if (dist2 <= radius*radius) {
           if (!IgnoreEqual || dist2 > 0.0)
               result.push_back(point);
       }

       // perpendicular distance to node boundary
       double perp = G::get(query, axis) - G::get(*point, axis);

       // change axis
       if (++axis >= K) {
           axis = 0;
       }

       // recurse to sub-trees if they are within radius
       if (sons[0] && perp <= +radius) {
           sons[0]->template range<IgnoreEqual>(query, radius, result, axis);
       }
       if (sons[1] && perp >= -radius) {
           sons[1]->template range<IgnoreEqual>(query, radius, result, axis);
       }
    }

    /**
     *  Returns all points that are within "radius" from "query".
     *  Returns points and their distance^2.
     */
    template<bool IgnoreEqual>
    void range(const T& query, double radius
               , std::vector<std::pair<T, double> > &result
               , unsigned int axis = 0) const
    {
        // return point if it is within radius
        T diff(G::diff(query, *point));
        double dist2 = inner_prod(diff, diff);
        if (dist2 <= radius*radius) {
            if (!IgnoreEqual || dist2 > 0.0)
                result.push_back(std::make_pair(*point, dist2));
        }

        // perpendicular distance to node boundary
        double perp = G::get(query, axis) - G::get(*point, axis);

        // change axis
        if (++axis >= K) {
            axis = 0;
        }

        // recurse to sub-trees if they are within radius
        if (sons[0] && perp <= +radius) {
            sons[0]->template range<IgnoreEqual>(query, radius, result, axis);
        }
        if (sons[1] && perp >= -radius) {
            sons[1]->template range<IgnoreEqual>(query, radius, result, axis);
        }
    }

    ~KdTreeNode() {
        delete sons[0];
        delete sons[1];
    }

protected:
    double inner_prod(const T &op1, const T &op2) const {
        double retval( 0.0 );
        for (unsigned int i = 0; i < K; i++) {
            retval += G::get(op1, i) * G::get(op2, i);
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

template <typename T>
struct GetCoordinate {
    typename T::value_type get(const T &value, unsigned int axis) const
    {
        return value(axis);
    }

    T diff(const T &op1, const T &op2) const {
        return op1 - op2;
    }
};

/**
 *  KdTree -- implementation of a k-d tree. The k-d tree is a generalization of
 *  the binary search tree which enables fast nearest neighbor searches in
 *  more dimensions (k > 1). See http://en.wikipedia.org/wiki/K-d_tree for
 *  more details.
 *
 *  typename T:     point in space
 *  unsinged int K: number of dimensions
 *  typename G:     dimension value getter (defaults to T::operator())
 *  typaneme C:     associated container type (defaults to std::vector)
 *
 *  NB: Associated container must outlive this tree (tree holds iterators to it)
 *      Associated container is not modified by any means.
 */
template<typename T, unsigned int K, typename G, typename C>
class KdTree {
private:
    typedef detail::KdTreeNode<T, K, G, C> node_type;
    typedef typename C::iterator mutable_iterator;

public:
    typedef typename C::const_iterator iterator;
    typedef std::size_t size_type;
    typedef std::pair<T, double> Neighbor;
    typedef std::vector<Neighbor> Neighbors;

    KdTree(const iterator &beg, const iterator &end)
        : begin_(beg), end_(end), size_(std::distance(beg, end))
    {
        if (!size_) { return; }

        typename node_type::Indirect indirect;
        indirect.reserve(std::distance(beg, end));
        for (iterator i(beg) ; i != end; ++i) {
            indirect.push_back(i);
        }

        root_.reset(new node_type(indirect.begin(), indirect.end()));
    }

    KdTree(const mutable_iterator &beg, const mutable_iterator &end
           , const IntrusiveKdTree &intrusive)
        : begin_(beg), end_(end), size_(std::distance(beg, end))
    {
        if (!size_) { return; }

        root_.reset(new node_type(beg, end, intrusive));
    }

    iterator begin() const { return begin_; }
    iterator end() const { return end_; }

    size_type size() const { return size_; }

    template<bool IgnoreEqual = false>
    iterator nearest(const T& query, double& dist2) const
    {
        if (!root_) {
            dist2 = .0;
            return end_;
        }
        return root_->template nearest<IgnoreEqual>(query, dist2);
    }

    template<bool IgnoreEqual = false>
    std::vector<T> range(const T& query, double radius) const {
        std::vector<T> result;
        if (!root_) { return result; }

        root_->template range<IgnoreEqual>(query, radius, result);
        return result;
    }

    template<bool IgnoreEqual = false>
    std::pair<iterator, double> nearest(const T& query) const
    {
        std::pair<iterator, double> res(end_, 0);
        if (!root_) { return res; }
        res.first = root_->template nearest<IgnoreEqual>(query, res.second);
        return res;
    }

    template<bool IgnoreEqual = false>
    void range(const T& query, double radius, std::vector<T>& result) const {
        if (!root_) { return; }
        return root_->template range<IgnoreEqual>(query, radius, result);
    }

    template<bool IgnoreEqual = false>
    void range(const T& query, double radius, std::vector<iterator>& result) const {
        if (!root_) { return; }
        return root_->template range<IgnoreEqual>(query, radius, result);
    }

    template<bool IgnoreEqual = false>
    void range(const T& query, double radius
               , std::vector<std::pair<T, double> > &result) const
    {
        if (!root_) { return; }
        return root_->template range<IgnoreEqual>(query, radius, result);
    }

private:
    boost::shared_ptr<node_type> root_;
    iterator begin_;
    iterator end_;
    size_type size_;
};

template<typename T, unsigned int K, typename G, typename C>
inline double samplingDelta(const KdTree<T, K, G, C> &tree
                            , double bulkThreshold)
{
    std::vector<double> distances;
    distances.reserve(std::distance(tree.begin(), tree.end()));
    for (const auto &point
             : boost::make_iterator_range(tree.begin(), tree.end()))
    {
        distances.push_back(tree.template nearest<true>(point).second);
    }

    std::sort(distances.begin(), distances.end());
    return std::sqrt(distances[distances.size() * bulkThreshold]);
}

} // namespace geometry

#endif

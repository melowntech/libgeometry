/**
 * Copyright (c) 2020 Melown Technologies SE
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
 * @file estimate-point-normals.hpp
 * @author Matyas Hollmann <matyas.hollmann@melowntech.com>
 *
 * Estimation of point normals in the point cloud.
 * Code based on window-mesh-legacy.cpp
 * 
 * Note: Uses Eigen3 library.
 *       
 */

#ifndef ESTIMATE_POINT_NORMALS_HPP_INCLUDED
#define ESTIMATE_POINT_NORMALS_HPP_INCLUDED

#include <Eigen/Dense>

#include "dbglog/dbglog.hpp"
#include "utility/openmp.hpp"

#include "kdtree.hpp"
#include "neighbors.hpp"

namespace geometry {
/**
 * Estimate the normal of a point (in K-dim space) based on the data matrix (N x K) 
 * which as rows has point itself and its (N - 1) closest neighbors. 
 * The estimation of the normal is based on the principal-component analysis
 * of the data, a SVD of the covariance matrix (K x K) is performed.
 **/
Eigen::VectorXd estimateNormal(const Eigen::MatrixXd& data);

/**
 *  Interface to access dimension values of a point, and to calculate the
 *  difference of two points.
 **/
template <typename T>
struct DefaultAccessor {
    typedef typename T::value_type value_type;
    typedef T difference_type;

    static inline typename T::value_type get(const T& pt, unsigned dim) {
        return pt(dim);
    }

    static inline void set(T& pt, unsigned dim, const typename T::value_type& val) {
        pt(dim) = val;
    }

    // used in Kdtree in calculation of the distance of two points
    static inline T diff(const T& op1, const T& op2) {
        return op1 - op2;
    }
};

/**
 *
 *  Estimate normals of all points in the input point pointcloud. 
 * 
 *  typename T:     point in a K-dimensional space
 *  unsigned K:     number of dimensions of the space
 *  typename A:     dimension value accessor (has to support get(),
 *                               set() and also diff() [used in Kdtree]
 */
template<typename T, unsigned K = 3, typename A = DefaultAccessor<T>>
std::vector<T> estimateNormals(const std::vector<T>& pointCloud,
                               unsigned nEstimatorPts = 40,
                               double radius = 0)
{
    static_assert(K >= 2,
                  "Estimation of point normals makes sense only in at least "
                  "2-dimensional space.");
    using Neighbor = typename KdTree<T, K, A>::Neighbor;
    using Neighbors = typename KdTree<T, K, A>::Neighbors;
    
    auto toEigen([](const Neighbor& n) -> Eigen::RowVectorXd {
        Eigen::RowVectorXd vec(K);
        for (unsigned t = 0; t < K; ++t) {
            vec(t) = A::get(n.first, t);
        }
        return vec;
    });
    
    auto fromEigen([](const Eigen::VectorXd& vec) -> T {
        T res;
        for (unsigned t = 0; t < K; ++t) {
            A::set(res, t, vec(t));
        }
        return res;
    });

    auto setZeros([](T& pt) -> void {
        for (unsigned t = 0; t < K; ++t) {
            A::set(pt, t, 0);
        }
    });

    const size_t nPoints(pointCloud.size());
    // prepare space for normals
    std::vector<T> normals(nPoints);
    
    LOG(info3) << "Building a kd-tree from the pointcloud of "
               << nPoints << " points.";
    KdTree<T, K, A> kdtree(pointCloud.begin(), pointCloud.end());

    /** per thread accumulative variables **/
    double searchRadiusTotal(0.0);
    size_t pointsProcessed(0);

    UTILITY_OMP(parallel for schedule(static) default(shared) 
               firstprivate(searchRadiusTotal, pointsProcessed))
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(nPoints); ++i)
    {
        const T& point(pointCloud[i]);
        T& normal(normals[i]);

        // Mode 1: use provided radius as search radius
        double searchRadius = radius;
        if (radius <= 0.0) {
            // Mode 2: search radius is variable and is based on the average
            // radius needed to reach the specified number of neighbors
            searchRadius = ((searchRadiusTotal && pointsProcessed)
                                ? (searchRadiusTotal / pointsProcessed)
                                : 1.0);
            ++pointsProcessed;
        }

        Neighbors neighbors;
        // find neighbors
        searchRadiusTotal += collectNeighbors(kdtree, point, neighbors
                                              , nEstimatorPts
                                              , searchRadius
                                              , (radius > 0));

        size_t nNeighs = neighbors.size();
        if (nNeighs < nEstimatorPts) {
            // Oops, normal stays zero -> we may deal with this later
            LOG(warn3) << "too few neighbors! (" << nNeighs << ")";
            setZeros(normal);
            continue;
        }
        
        Eigen::MatrixXd samples(nNeighs, K);
        for (size_t j = 0; j < nNeighs; ++j)
        {
            samples.row(j) = toEigen(neighbors[j]);
        }

        // calculate normal
        normal = fromEigen(estimateNormal(samples));
    }
    return normals;
}

/**
 * @brief Orients point cloud normals to outward directions.
 *
 * @param pointCloud  Input points cloud
 * @param normals     Estimated normals (with undetermined orientation) for each point
 * @param pointRadius Radius of a patch associated with each point. Must be large enough
 *                    for the point patches to overlap and cover the surface.
 */
void reorientNormals(const std::vector<math::Point3>& pointCloud
                     , std::vector<math::Point3>& normals
                     , double pointRadius);

} // geometry

#endif // ESTIMATE_POINT_NORMALS_HPP_INCLUDED

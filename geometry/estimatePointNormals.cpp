#include <Eigen/Dense>

#include "dbglog/dbglog.hpp"
#include "utility/openmp.hpp"

#include "./kdtree.hpp"

namespace geometry {
/**
 * Estimate the normal of a point (in K-dim space) based on the data matrix (N x K) 
 * which as rows has point itself and its (N - 1) closest neighbors. 
 * The estimation of the normal is based on the principal-component analysis
 * of the data, a SVD of the covariance matrix (K x K) is performed.
 **/
Eigen::VectorXd estimateNormal(const Eigen::MatrixXd& data)
{
    using namespace Eigen;
    // calculate centroid (1 x K) of all (N) data points
    RowVectorXd centroid = data.colwise().mean();
    // mean-center sample matrix
    MatrixXd centered = data.rowwise() - centroid;
    // calculate the covariance matrix (K x K), for PCA we may ommit
    // scaling by the factor 1/(N - 1)
    MatrixXd covariance = centered.transpose() * centered;

    // calculate SVD (ordered by the magnitude of SV)
    JacobiSVD<MatrixXd> svd(covariance, ComputeFullU);
    // the normal is the last column of U, singular vector of the covariance
    // matrix corresponding to the smallest singular value 
    // (i.e., the direction of the lowest variability in the data)
    VectorXd sgVec(svd.matrixU().col(data.cols() - 1));

    return sgVec.normalized();
}

/**
 *  Getter/setter of 
 * 
 **/
template <typename T>
struct DefaultAccessor {
    static inline typename T::value_type get(const T& pt, unsigned dim) {
        return pt(dim);
    }

    static inline void set(T& pt, unsigned dim, typename T::value_type val) {
        pt(dim) = val;
    }

    // used in Kdtree in calculation of the distance of two points
    static inline T diff(const T& op1, const T& op2) {
        return op1 - op2;
    }
};

/**
 *  typename T:     point in K-dimensional space
 *  unsinged K:     number of dimensions of the space
 *  typename A:     dimension value accessor (has to support get(),
 *                               set() and diff() (see kdtree.hpp))
 */
template<typename T, unsigned K = 3, typename A = DefaultAccessor<T>>
std::vector<T> estimateNormals(const std::vector<T>& pointCloud,
                               unsigned nEstimatorPts = 20,
                               double radius = 0)
{
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

    const size_t nPoints(pointCloud.size());
    // prepare space for normals
    std::vector<T> normals(nPoints, T());
    
    LOG(info3) << "Building a kd-tree from the pointcloud of "
               << nPoints << " points.";
    KdTree<T, K, A> kdtree(pointCloud.begin(), pointCloud.end());

    /** per thread variables*/
    double searchRadiusTotal(0.0);
    size_t pointsProcessed(0);

   UTILITY_OMP(parallel for schedule(static) default(shared) 
               firstprivate(searchRadiusTotal, pointsProcessed))
    for (size_t i = 0; i < nPoints; ++i)
    {
        const T& point(pointCloud[i]);
        T& normal(normals[i]);

        // Mode 1: use provided radius as search radius
        double searchRadius = radius;
        if (radius <= 0.0) {
            // Mode 2: search radius is variable and is based on the average
            // radius needed to reach the specified number of neighbors
            searchRadius
                = (pointsProcessed ? (searchRadiusTotal / pointsProcessed)
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
            // Oops, normal stays zero -> we ignore point later
            LOG(warn3) << "too few neighbors! (" << nNeighs << ")";
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
}

} // geometry
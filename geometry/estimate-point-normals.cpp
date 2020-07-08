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
#include "estimate-point-normals.hpp"
#include "bvh.hpp"
#include "math/math.hpp"
#include "math/geometry.hpp"
#include "utility/expect.hpp"

#include <boost/iterator/counting_iterator.hpp>

namespace geometry {

namespace ublas = math::ublas;

Eigen::VectorXd estimateNormal(const Eigen::MatrixXd& data
                               , const bool confidence)
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
    // (i.e., the direction of the least variability in the data)
    VectorXd sgVec(svd.matrixU().col(data.cols() - 1));

    double weight = 1.;
    if (confidence) {
        double norm = svd.singularValues().norm();
        weight = 1. - svd.singularValues()[data.cols() - 1] / norm;
        weight = (exp(8. * weight) - 1.) / (exp(8.) - 1.);
    }
    // normalization shouldn't be needed as U should be a unitary matrix.
    return sgVec * weight;
}

namespace {

    class BvhDisk : public BvhPrimitive {
        math::Point3 center_;
        math::Point3 normal_;
        double radius_;

    public:
        BvhDisk() = default;

        BvhDisk(const math::Point3& center
                , const math::Point3& normal
                , const double radius
                , const std::uint32_t index)
            : center_(center)
            , normal_(normal)
            , radius_(radius) {
            assert(radius_ > 0.);
            userData = index;
        }

        bool getIntersection(const Ray& ray, IntersectionInfo& intersection) const {
            // ray-plane intersection
            const double denom = inner_prod(normal_, ray.direction());
            if (std::fabs(denom) < 1.e-6) {
                return false;
            }
            const math::Point3 diff = center_ - ray.origin();
            const double t = inner_prod(diff, normal_) / denom;
            if (t > 0.) {
                // find the distance of the intersection from the disk center
                const math::Point3 dp = ray.origin() + ray.direction() * t - center_;
                const double distSqr = inner_prod(dp, dp);
                if (distSqr < math::sqr(radius_)) {
                    intersection.object = this;
                    intersection.t = t;
                    return true;
                }
            }

            return false;
        }

        math::Extents3 getBBox() const {
            ublas::scalar_vector<double> size(3, radius_);
            return math::Extents3(center_ - size, center_ + size);
        }

        math::Point3 getCenter() const {
            return center_;
        }
    };

} // namespace

static constexpr float DX = 0.25;
static constexpr float DZ = 1.;
static math::Points3 ALL_DIRS = {
    math::Point3(0., 0., DZ),
    math::Point3(0., DX, DZ),
    math::Point3(0., -DX, DZ),
    math::Point3(DX, 0., DZ),
    math::Point3(-DX, 0., DZ),
    math::Point3(DX, DX, DZ),
    math::Point3(DX, -DX, DZ),
    math::Point3(-DX, -DX, DZ),
    math::Point3(-DX, DX, DZ),
};

void reorientNormals(const std::vector<math::Point3>& pc
                     , std::vector<math::Point3>& normals
                     , const double pointRadius)
{
    std::vector<BvhDisk> disks(pc.size());
    UTILITY_OMP(parallel for)
    for (std::uint32_t i = 0; i < pc.size(); ++i) {
        disks[i] = BvhDisk(pc[i], normals[i], pointRadius, i);
    }

    LOG(info2) << "Building BVH for " << disks.size() << " points";
    geometry::Bvh<BvhDisk> bvh(20);
    bvh.build(std::move(disks));

    // find order in z-direction
    std::vector<std::uint32_t> orderZ(pc.size());
    std::copy(boost::counting_iterator<std::uint32_t>(0)
              , boost::counting_iterator<std::uint32_t>(pc.size())
              , orderZ.begin());
    std::sort(orderZ.begin(), orderZ.end()
              , [&](std::uint32_t i1, std::uint32_t i2) { return pc[i1](2) > pc[i2](2); });

    // step 1 - determine normal orientation based on normals of already determined points,
    // assuming the topmost points (roofs) have z>0 orientation.
    LOG(info2) << "Estimating normal orientations";
    UTILITY_OMP(parallel for)
    for (std::uint32_t rankZ = 0; rankZ < orderZ.size(); ++rankZ) {
        const std::uint32_t i = orderZ[rankZ];
        int doFlip = 0;
        for (math::Point3 dir : ALL_DIRS) {
            dir = math::normalize(dir);
            // whether the ray is outward or inward with respect to current normal orientation
            const int outward = math::sgn(inner_prod(dir, normals[i]));
            const Ray ray(pc[i] + pointRadius * dir, dir);
            IntersectionInfo is;
            if (bvh.getFirstIntersection(ray, is)) {
                const std::uint32_t j = is.object->userData;
                // if this is an outward ray, flip the normal if the intersected point has
                // the same orientation as the ray (i.e. likely a backface);
                // for an inward ray, flip the normal if the intersected point is NOT a backface.
                doFlip += int(inner_prod(normals[j], dir) > 0.) * outward;
            } else {
                // no occlusion -> normal correct if this is an outward ray
                doFlip -= outward;
            }
        }

        if (doFlip > 0) {
            // more votes for flipping the normal
            normals[i] *= -1;
        }
    }

    // step 2 - flip normals which have opposite orientation than their neighbors
    LOG(info2) << "Flipping outliers";
    KdTree<math::Point3, 3> tree(pc.begin(), pc.end());
    std::vector<uint8_t> doFlip(pc.size(), false);
    std::vector<math::Points3::const_iterator> neighs;
    UTILITY_OMP(parallel for private(neighs))
    for (std::uint32_t i = 0; i < pc.size(); ++i) {
        neighs.clear();
        tree.range(pc[i], 2 * pointRadius, neighs);

        double count = 0.;
        for (auto& n : neighs) {
            const std::uint32_t j = std::uint32_t(n - pc.begin());
            count += inner_prod(normals[i], normals[j]);
        }
        if (count < 0.) {
            doFlip[i] = true;
        }
    }
    UTILITY_OMP(parallel for)
    for (std::uint32_t i = 0; i < pc.size(); ++i) {
        if (doFlip[i]) {
            normals[i] *= -1;
        }
    }
}

} // namespace geometry

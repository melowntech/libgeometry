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

#include <boost/iterator/counting_iterator.hpp>

namespace geometry {

namespace ublas = math::ublas;

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
    // (i.e., the direction of the least variability in the data)
    VectorXd sgVec(svd.matrixU().col(data.cols() - 1));

    // normalization shouldn't be needed as U should be a unitary matrix.
    return sgVec;//.normalized();
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
                , const uint index)
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
    for (uint i = 0; i < pc.size(); ++i) {
        disks[i] = BvhDisk(pc[i], normals[i], pointRadius, i);
    }

    LOG(info2) << "Building BVH for " << disks.size() << " points";
    geometry::Bvh<BvhDisk> bvh(20);
    bvh.build(std::move(disks));

    // find order in z-direction
    std::vector<uint> orderZ(pc.size());
    std::copy(boost::counting_iterator<uint>(0)
              , boost::counting_iterator<uint>(pc.size())
              , orderZ.begin());
    std::sort(orderZ.begin(), orderZ.end()
              , [&](uint i1, uint i2) { return pc[i1](2) > pc[i2](2); });

    LOG(info2) << "Estimating normal orientations";
    UTILITY_OMP(parallel for)
    for (uint rankZ = 0; rankZ < orderZ.size(); ++rankZ) {
        const uint i = orderZ[rankZ];
        int doFlip = 0;
        for (math::Point3 dir : ALL_DIRS) {
            dir = math::normalize(dir);
            // whether the ray is outward or inward with respect to current normal orientation
            const int outward = math::sgn(inner_prod(dir, normals[i]));
            const Ray ray(pc[i] + pointRadius * dir, dir);
            IntersectionInfo is;
            if (bvh.getFirstIntersection(ray, is)) {
                const uint j = is.object->userData;
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
}

} // namespace geometry

/**
 * Copyright (c) 2019 Melown Technologies SE
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
 * @file bvh.hpp
 * @author Pavel Sevecek <pavel.sevecek@melown.com>
 *
 * Bounding volume hierarchy
 */

#ifndef bvh_hpp_included_
#define bvh_hpp_included_

#include "math/geometry_core.hpp"

namespace geometry {

class Ray {
    friend bool intersectBox(const math::Extents3& box, const Ray& ray, double& t_min, double& t_max);

    math::Point3 orig_;
    math::Point3 dir_;
    math::Point3 invDir_;
    std::array<int, 3> signs_;

public:
    Ray() = default;

    Ray(const math::Point3& origin, const math::Point3& dir)
        : orig_(origin)
        , dir_(dir) {
        for (int i = 0; i < 3; ++i) {
            invDir_[i] = (dir[i] == 0.) ? INFINITY : 1. / dir[i];
            signs_[i] = int(invDir_[i] < 0.);
        }
    }

    const math::Point3& origin() const {
        return orig_;
    }

    const math::Point3& direction() const {
        return dir_;
    }
};

inline bool intersectBox(const math::Extents3& box, const Ray& ray, double& t_min, double& t_max) {
    std::array<math::Point3, 2> b{ box.ll, box.ur };
    double tmin = (b[ray.signs_[0]](0) - ray.orig_(0)) * ray.invDir_(0);
    double tmax = (b[1 - ray.signs_[0]](0) - ray.orig_(0)) * ray.invDir_(0);
    assert(!std::isnan(tmin) && !std::isnan(tmax)); // they may be inf though
    const double tymin = (b[ray.signs_[1]](1) - ray.orig_(1)) * ray.invDir_(1);
    const double tymax = (b[1 - ray.signs_[1]](1) - ray.orig_(1)) * ray.invDir_(1);
    assert(!std::isnan(tymin) && !std::isnan(tymax));

    if ((tmin > tymax) || (tymin > tmax)) {
        return false;
    }
    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);

    const double tzmin = (b[ray.signs_[2]](2) - ray.orig_(2)) * ray.invDir_(2);
    const double tzmax = (b[1 - ray.signs_[2]](2) - ray.orig_(2)) * ray.invDir_(2);
    assert(!std::isnan(tzmin) && !std::isnan(tzmax));

    if ((tmin > tzmax) || (tzmin > tmax)) {
        return false;
    }
    tmin = std::max(tmin, tzmin);
    tmax = std::min(tmax, tzmax);

    t_min = tmin;
    t_max = tmax;

    return true;
}


struct BvhPrimitive {
    /** Generic user data, can be used to store additional information to the primitives. */
    uint userData = uint(-1);
};

/**
 * \brief Holds intormation about intersection.
 */
struct IntersectionInfo {
    /** Distance of the hit in units of ray.direction(). */
    double t;

    /** Object hit by the ray, or nullptr if nothing has been hit. */
    const BvhPrimitive* object = nullptr;

    math::Point3 hit(const Ray& ray) const {
        return ray.origin() + ray.direction() * t;
    }

    /** Sort by intersection distance. */
    bool operator<(const IntersectionInfo& other) const {
        return t < other.t;
    }
};

/**
 * @brief Bounding volume hierarchy
 * @details TBvhObject must derive from \ref BvhPrimitive and implement the following interface:
 * @code
 * struct BvhObject : public BvhPrimitive {
 *   bool getIntersection(const Ray& ray, IntersectionInfo& intersection) const;
 *   math::Extents3 getBBox() const;
 *   math::Point3 getCenter() const;
 * };
 * @endcode
 */
template <typename TBvhObject>
class Bvh : public boost::noncopyable {
private:
    const uint leafSize_;
    uint nodeCnt_ = 0;
    uint leafCnt_ = 0;

    std::vector<TBvhObject> objects_;

    struct BvhNode {
        math::Extents3 box;
        uint start;
        uint primCnt;
        uint rightOffset;
    };

    std::vector<BvhNode> nodes_;

    struct BvhTraversal {
        uint idx;
        double t_min;
    };

public:
    explicit Bvh(const uint leafSize = 4)
        : leafSize_(leafSize) {}

    /// \brief Contructs the BVH from given set of objects.
    ///
    /// This erased previously stored objects.
    void build(std::vector<TBvhObject>&& objects) {
        objects_ = std::move(objects);
        assert(!objects_.empty());
        nodeCnt_ = 0;
        leafCnt_ = 0;

        struct BvhBuildEntry {
            uint parent;
            uint start;
            uint end;
        };

        std::array<BvhBuildEntry, 128> stack;
        uint stackIdx = 0;
        constexpr uint NO_PARENT_FLAG = uint(-1);
        constexpr uint UNTOUCHED_FLAG = uint(-1);

        // Push the root
        stack[stackIdx].start = 0;
        stack[stackIdx].end = objects_.size();
        stack[stackIdx].parent = NO_PARENT_FLAG;
        stackIdx++;

        BvhNode node;
        std::vector<BvhNode> buildNodes;
        buildNodes.reserve(2 * objects_.size());

        while (stackIdx > 0) {
            BvhBuildEntry& nodeEntry = stack[--stackIdx];
            const uint start = nodeEntry.start;
            const uint end = nodeEntry.end;
            const uint primCnt = end - start;

            nodeCnt_++;
            node.start = start;
            node.primCnt = primCnt;
            node.rightOffset = UNTOUCHED_FLAG;

            math::Extents3 bbox = objects_[start].getBBox();
            const math::Point3 center = objects_[start].getCenter();
            math::Extents3 boxCenter(center, center);
            for (uint i = start + 1; i < end; ++i) {
                math::update(bbox, objects_[i].getBBox());
                math::update(boxCenter, objects_[i].getCenter());
            }
            node.box = bbox;

            if (primCnt <= leafSize_) {
                node.rightOffset = 0;
                leafCnt_++;
            }
            buildNodes.push_back(node);

            if (nodeEntry.parent != NO_PARENT_FLAG) {
                buildNodes[nodeEntry.parent].rightOffset--;

                if (buildNodes[nodeEntry.parent].rightOffset == UNTOUCHED_FLAG - 2) {
                    buildNodes[nodeEntry.parent].rightOffset = nodeCnt_ - 1 - nodeEntry.parent;
                }
            }

            if (node.rightOffset == 0) {
                continue;
            }

            const uint splitDim = argMax(boxCenter.ur - boxCenter.ll);
            const double split = 0.5 * (boxCenter.ll(splitDim) + boxCenter.ur(splitDim));

            uint mid = start;
            for (uint i = start; i < end; ++i) {
                if (objects_[i].getCenter()[splitDim] < split) {
                    std::swap(objects_[i], objects_[mid]);
                    ++mid;
                }
            }

            if (mid == start || mid == end) {
                mid = start + (end - start) / 2;
            }

            if (stackIdx >= stack.size() - 2) {
                objects_.clear();
                nodeCnt_ = leafCnt_ = 0;
                throw std::runtime_error("BVH build stack overflow, try increasing the leaf size");
            }

            stack[stackIdx++] = { nodeCnt_ - 1, mid, end };
            stack[stackIdx++] = { nodeCnt_ - 1, start, mid };
        }

        assert(buildNodes.size() == nodeCnt_);
        nodes_ = std::move(buildNodes);
    }

    /// \brief Finds the closest intersection of the ray.
    ///
    /// Returns true if an intersection has been found.
    bool getFirstIntersection(const Ray& ray, IntersectionInfo& intersection) const {
        intersection.t = INFINITY;
        intersection.object = nullptr;

        this->getIntersections(ray, [&intersection](IntersectionInfo& current) {
            if (current.t < intersection.t) {
                intersection = current;
            }
            return true;
        });
        return intersection.object != nullptr;
    }

    /// \brief Returns all intersections of the ray.
    void getAllIntersections(const Ray& ray, std::set<IntersectionInfo>& intersections) const {
        intersections.clear();
        this->getIntersections(ray, [&intersections](IntersectionInfo& current) { //
            if (current.t > 0.) {
                intersections.insert(current);
            }
            return true;
        });
    }

    /// \brief Returns true if the ray is occluded by some geometry
    bool isOccluded(const Ray& ray) const {
        bool occluded = false;
        this->getIntersections(ray, [&occluded](IntersectionInfo&) {
            occluded = true;
            return false; // do not continue with traversal
        });
        return occluded;
    }

    /// \brief Returns the bounding box of all objects in BVH.
    math::Extents3 getBoundingBox() const {
        return nodes_[0].box;
    }

private:
    template <typename TAddIntersection>
    void getIntersections(const Ray& ray, const TAddIntersection& addIntersection) const {
        std::array<double, 4> boxHits;
        uint closer;
        uint other;

        std::array<BvhTraversal, 64> stack;
        int stackIdx = 0;

        stack[stackIdx].idx = 0;
        stack[stackIdx].t_min = 0.;

        while (stackIdx >= 0) {
            const uint idx = stack[stackIdx].idx;
            stackIdx--;
            const BvhNode& node = nodes_[idx];

            if (node.rightOffset == 0) {
                // leaf
                for (uint primIdx = 0; primIdx < node.primCnt; ++primIdx) {
                    IntersectionInfo current;

                    const TBvhObject& obj = objects_[node.start + primIdx];
                    const bool hit = obj.getIntersection(ray, current);

                    if (hit) {
                        if (!addIntersection(current)) {
                            // bailout
                            return;
                        }
                    }
                }
            } else {
                // inner node
                const bool hitLeft = intersectBox
                        (nodes_[idx + 1].box, ray, boxHits[0], boxHits[1]);
                const bool hitRight = intersectBox
                        (nodes_[idx + node.rightOffset].box, ray, boxHits[2], boxHits[3]);

                if (stackIdx >= int(stack.size()) - 2) {
                    throw std::runtime_error("BVH traversal stack overflow, try increasing the leaf size");
                }

                if (hitLeft && hitRight) {
                    closer = idx + 1;
                    other = idx + node.rightOffset;

                    if (boxHits[2] < boxHits[0]) {
                        std::swap(boxHits[0], boxHits[2]);
                        std::swap(boxHits[1], boxHits[3]);
                        std::swap(closer, other);
                    }

                    stack[++stackIdx] = BvhTraversal{ other, boxHits[2] };
                    stack[++stackIdx] = BvhTraversal{ closer, boxHits[0] };
                } else if (hitLeft) {
                    stack[++stackIdx] = BvhTraversal{ idx + 1, boxHits[0] };
                } else if (hitRight) {
                    stack[++stackIdx] = BvhTraversal{ idx + node.rightOffset, boxHits[2] };
                }
            }
        }
    }

    inline int argMax(const math::Point3& v) {
        int maxIdx = 0;
        if (v(1) > v(maxIdx)) {
            maxIdx = 1;
        }
        if (v(2) > v(maxIdx)) {
            maxIdx = 2;
        }
        return maxIdx;
    }
};

} // namespace geometry

#endif

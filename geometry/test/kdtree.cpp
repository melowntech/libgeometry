/**
 * Copyright (c) 2017 Melown Technologies SE
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
#include <boost/test/unit_test.hpp>

#include <math/math_all.hpp>
#include <geometry/kdtree.hpp>

BOOST_AUTO_TEST_CASE(geometry_kdtree_nearest)
{
    const int N = 1000;
    BOOST_TEST_MESSAGE("* Testing kd-tree on " << N << " points");

    // generate random points in the unit cube
    srand(0);
    std::vector<math::Point3> points;
    for (int i = 0; i < N; i++)
    {
        points.emplace_back((double) rand() / RAND_MAX,
                            (double) rand() / RAND_MAX,
                            (double) rand() / RAND_MAX);
    }

    // build kd-tree
    geometry::KdTree<math::Point3> kdtree(points.begin(), points.end());

    // find nearest neighbor for each point
    for (int i = 0; i < N; i++)
    {
        double dist2;
        const math::Point3& neigh = *kdtree.nearest<true>(points[i], dist2);

        // check with brute force calculation
        double bfDist2 = INFINITY;
        const math::Point3* bfNeigh = 0;
        for (int j = 0; j < N; j++)
        {
            if (j == i) continue;
            math::Point3 diff(points[i] - points[j]);
            double d2 = inner_prod(diff, diff);
            if (d2 < bfDist2)
            {
                bfDist2 = d2;
                bfNeigh = &points[j];
            }
        }

        // test
        BOOST_REQUIRE((dist2 == bfDist2) && (neigh == *bfNeigh));
    }
}



    // TODO: test kdtree.nearest<false>
    // TODO: test const iterator version

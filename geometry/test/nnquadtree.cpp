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

#include "math/math_all.hpp"
#include "geometry/nnquadtree.hpp"

BOOST_AUTO_TEST_CASE(geometry_nnquadtree)
{
    const int N = 10000;
    BOOST_TEST_MESSAGE("* Testing NNQuadTree on " << N << " points");

    // generate random points in the unit square
    math::Points2 points;
    srand(0);
    for (int i = 0; i < N; i++)
    {
        points.emplace_back((double) rand() / RAND_MAX,
                            (double) rand() / RAND_MAX);
    }

    // insert the points into a quadtree
    geometry::NNQuadTree<math::Point2, int, 16> qtree({0,0}, {1,1});
    for (int i = 0; i < N; i++)
    {
        qtree.insert(points[i], i);
    }

    // test NN search, compare with brute force algorithm
    for (int i = 0; i < N; i++)
    {
        math::Point2 testme((double) rand() / RAND_MAX,
                            (double) rand() / RAND_MAX);

        int id = -1;
        double minDist = std::numeric_limits<double>::max();
        for (int j = 0; j < N; j++)
        {
            double d = norm_2(points[j] - testme);
            if (d < minDist) {
                minDist = d;
                id = j;
            }
        }

        // test
        BOOST_REQUIRE(id == qtree.nearest(testme));
    }
}

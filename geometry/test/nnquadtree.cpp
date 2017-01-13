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

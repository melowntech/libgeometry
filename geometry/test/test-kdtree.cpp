#include <math/math_all.hpp>
#include <geometry/kdtree.hpp>



int main(int argc, char* argv[])
{
    const int N = 1000;
    std::cout << "Testing kd-tree on " << N << " points\n";

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
        if (dist2 != bfDist2 || neigh != *bfNeigh)
        {
            std::cout << "TEST FAILED (on point " << i << " of " << N << ")\n";
            return EXIT_SUCCESS;
        }
    }


    // TODO: test kdtree.nearest<false>
    // TODO: test const iterator version

    return EXIT_SUCCESS;
}


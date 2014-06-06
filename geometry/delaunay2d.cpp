/**
 * This is workaround for CGAL Assertions currently ruining debug mode because
 * CGAL expects flags which are not currently set
 */
#ifndef NDEBUG
#  define NDEBUG
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include "delaunay2d.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel            Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                    Delaunay;
typedef Kernel::Point_2                                                Point;

namespace geometry {

std::vector<DTriangle> delaunayTriangulation2d(const math::Points2 points)
{
    std::vector<std::pair<Point, unsigned> > cgalPoints;
    {
        unsigned index = 0;
        for (const auto &pt : points) {
            cgalPoints.emplace_back(Point(pt(0), pt(1)), index++);
        }
    }

    Delaunay triangulation;
    triangulation.insert(cgalPoints.begin(), cgalPoints.end());

    std::vector<DTriangle> result;
    result.reserve(triangulation.number_of_faces());

    for (auto fit = triangulation.finite_faces_begin();
              fit != triangulation.finite_faces_end(); ++fit)
    {
        DTriangle indices;
        for (int i = 0; i < 3; i++) {
            indices[i] = fit->vertex(i)->info();
        }
        result.push_back(indices);
    }

    return result;
}

} // namespace geometry

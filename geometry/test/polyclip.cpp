#include <cstdlib>

#include <boost/test/unit_test.hpp>

#include "geometry/polyclip.hpp"

#include "dbglog/dbglog.hpp"

namespace {
    math::Viewport2 viewport(200, 100, 50, 55);

    std::vector<math::Points2i> polygons
    ({
        {{50, 0}, {280, 70}, {170, 190}, {5, 60}}
        , {{50, 0}, {280, 70}, {170, 190}, {5, 105}}
        , {{50, 0}, {280, 70}, {170, 190}, {25, 5}}
        , {{25, 5}, {170, 190}, {280, 70}, {50, 0}}
        , {{50, 0}, {280, 70}, {170, 150}, {5, 60}}
        , {{5, 5}, {225, 5}, {225, 115}, {5, 115}}
        , {{-100, 60}, {120, -100}, {320, 60}, {120, 220}}
        , {{25, 5}, {280, 70}, {170, 190}, {50, 0}}

        , {{150, 0}, {280, 70}, {170, 200}, {5, 60}}

        , {{150, 0}, {0, 150}, {0, 0}}
        , {{150, 0}, {297, 150}, {297, 0}}

        , {{150, 210}, {0, 60}, {0, 210}}

        , {{150, 210}, {0, 150}, {0, 210}}
    });

    std::vector<math::Points2i> clippedPolygons
    ({
        {{50, 55}, {230, 55}, {250, 60}, {250, 102}
            , {202, 155}, {125, 155}, {50, 95}}
        , {{50, 55}, {230, 55}, {250, 60}, {250, 102}
            , {202, 155}, {102, 155}, {50, 128}}
        , {{230, 55}, {250, 60}, {250, 102}, {202, 155}
            , {142, 155}, {64, 55}}
        , {{64, 55}, {142, 155}, {202, 155}, {250, 102}
            , {250, 60}, {230, 55}}
        , {{50, 55}, {230, 55}, {250, 60}, {250, 91}, {170, 150}
            , {50, 84}}
        , {{50, 55}, {225, 55}, {225, 115}, {50, 115}}
        , {{50, 155}, {50, 55}, {250, 55}, {250, 116}, {201, 155}}
        , {{221, 55}, {250, 62}, {250, 102}, {202, 155}
            , {147, 155}, {84, 55}}
        , {{50, 55}, {250, 55}, {250, 105}, {208, 155}, {116, 155}, {50, 98}}
        , {{50, 55}, {95, 55}, {50, 100}}
        , {{250, 55}, {203, 55}, {250, 102}}
        , {{50, 155}, {95, 155}, {50, 110}}
        , {}
    });

} // namespace

template <typename T>
void printPolygon(const std::vector<math::Point2_<T> > &points)
{
    std::cout << "(";
    for (const auto &p : points) {
        std::cout << "(" << p(0) << ", " << p(1) << "), ";
    }
    std::cout << ")" << std::endl;
}

void printViewport(const math::Viewport2 &v)
{
    std::cout << "viewport = (" << v.x << ", " << v.y
              << ", " << v.width << ", " << v.height << ")" << std::endl;
}

BOOST_AUTO_TEST_CASE(geometry_polyclip_1)
{
    BOOST_TEST_MESSAGE("* Testing polygon clipping.");

    std::vector<math::Points2i> clipped;
    for (const auto &polygon : polygons) {
        clipped.push_back(geometry::clip(viewport, polygon));
    }

    int idx(0);
    auto iclippedPolygons(clippedPolygons.cbegin());
    for (const auto &c : clipped) {
        bool result(std::equal
                    (c.begin(), c.end(), iclippedPolygons->begin()));
        BOOST_REQUIRE_MESSAGE(result, "Result differs: " << idx);
        ++iclippedPolygons;
        ++idx;
    }

    if (std::getenv("DEBUG_POLYGON_CLIPPING")) {
        printViewport(viewport);

        std::cout << "data = (" << std::endl;
        auto iclipped(clipped.cbegin());

        for (const auto &polygon : polygons) {
            std::cout << "{ \"polygon\": ";
            printPolygon(polygon);
            std::cout << ", \"clipped\": ";
            printPolygon(*iclipped++);
            std::cout << "}, " << std::endl;
        }
        std::cout << ")" << std::endl;
    }
}

define_module(LIBRARY geometry DEPENDS math dbglog OPENMESH)

set(geometry_SOURCES
	pointcloud.cpp
	volume.cpp
    parse-obj.cpp

    polygon.hpp
    polyclip.cpp

    mesh.hpp
    mesh.cpp
)

add_library(geometry STATIC ${geometry_SOURCES})
buildsys_library(geometry)
target_link_libraries(geometry ${MODULE_LIBRARIES})

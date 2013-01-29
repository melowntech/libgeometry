define_module(LIBRARY geometry DEPENDS dbglog)

set(geometry_SOURCES
	pointcloud.cpp
	volume.cpp
    parse-obj.cpp

    polygon.hpp
    polyclip.cpp

    mesh.hpp
)

add_library(geometry STATIC ${geometry_SOURCES})
set_target_properties(geometry PROPERTIES ARCHIVE_OUTPUT_DIRECTORY
  ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(geometry ${MODULE_LIBRARIES})

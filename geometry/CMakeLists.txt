define_module(LIBRARY geometry DEPENDS math utility dbglog OPENMESH)

if (OPENMESH_FOUND)
  message(STATUS "geometry: compiling in OpenMesh support")
  append_parent_directory_list_property(COMPILE_DEFINITIONS
    "GEOMETRY_HAS_OPENMESH=1")

  set(geometry_OPENMESH_SOURCES
    meshop-openmesh.cpp
    )
else()
  message(STATUS "geometry: compiling without OpenMesh support")
endif()

set(geometry_SOURCES
	pointcloud.cpp
	volume.cpp
    parse-obj.cpp

    polygon.hpp
    polyclip.cpp

    mesh.hpp mesh.cpp
    meshop.hpp meshop.cpp

    ${geometry_OPENMESH_SOURCES}
)

add_library(geometry STATIC ${geometry_SOURCES})
buildsys_library(geometry)
target_link_libraries(geometry ${MODULE_LIBRARIES})

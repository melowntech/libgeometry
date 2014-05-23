set(geometry_EXTRA_DEPENDS)

if (OPENMESH_FOUND)
  message(STATUS "geometry: compiling in OpenMesh support")
  append_parent_directory_list_property(COMPILE_DEFINITIONS
    "GEOMETRY_HAS_OPENMESH=1")

  set(geometry_OPENMESH_SOURCES
    meshop-openmesh.cpp
    )

  list(APPEND geometry_EXTRA_DEPENDS OPENMESH)
else()
  message(STATUS "geometry: compiling without OpenMesh support")
endif()

if (CGAL_FOUND)
  message(STATUS "geometry: compiling in CGAL support")
  append_parent_directory_list_property(COMPILE_DEFINITIONS
    "GEOMETRY_HAS_CGAL=1")

  set(geometry_CGAL_SOURCES
    delaunay2d.cpp delaunay2d.hpp
    )

  list(APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} ${GMP_LIBRARIES}) # TODO: hide
  list(APPEND geometry_EXTRA_DEPENDS CGAL)
else()
  message(STATUS "geometry: compiling without CGAL support")
endif()


define_module(LIBRARY geometry DEPENDS math utility dbglog
  ${geometry_EXTRA_DEPENDS})

set(geometry_SOURCES
	pointcloud.cpp
	volume.cpp
    parse-obj.cpp

    polygon.hpp
    polyclip.cpp

    mesh.hpp mesh.cpp
    meshop.hpp meshop.cpp

    ${geometry_OPENMESH_SOURCES}
    ${geometry_CGAL_SOURCES}
)

add_library(geometry STATIC ${geometry_SOURCES})
buildsys_library(geometry)
target_link_libraries(geometry ${MODULE_LIBRARIES})

# bump version here
set(geometry_VERSION 1.7)

set(geometry_EXTRA_DEPENDS)
set(geometry_EXTRA_DEFINITIONS)

if (OPENMESH_FOUND)
  message(STATUS "geometry: compiling in OpenMesh support")

  # add definitions
  list(APPEND geometry_EXTRA_DEFINITIONS GEOMETRY_HAS_OPENMESH=1)

  set(geometry_OPENMESH_SOURCES
    meshop-openmesh.cpp
    )

  # handle() is deprecated -> no warning
  set_source_files_properties(meshop-openmesh.cpp PROPERTIES
    COMPILE_FLAGS -Wno-deprecated-declarations)

  list(APPEND geometry_EXTRA_DEPENDS OPENMESH)
else()
  message(STATUS "geometry: compiling without OpenMesh support")
endif()

if (CGAL_FOUND)
  message(STATUS "geometry: compiling in CGAL support")

  # add definitions
  list(APPEND geometry_EXTRA_DEFINITIONS GEOMETRY_HAS_CGAL=1)

  set(geometry_CGAL_SOURCES
    delaunay2d.cpp delaunay2d.hpp
    )

  list(APPEND geometry_EXTRA_DEPENDS CGAL)
else()
  message(STATUS "geometry: compiling without CGAL support")
endif()

if (OpenCV_FOUND)
  message(STATUS "geometry: compiling in OpenCV support")

  # add definitions
  list(APPEND geometry_EXTRA_DEFINITIONS GEOMETRY_HAS_OPENCV=1)

  set(geometry_OPENCV_SOURCES
    faceclip.hpp faceclip.cpp
    pointindex.hpp
    )

  list(APPEND geometry_EXTRA_DEPENDS OpenCV)
else()
  message(STATUS "geometry: compiling without OpenCV support")
endif()


define_module(LIBRARY geometry=${geometry_VERSION}
  DEPENDS math>=1.2 utility>=1.11 dbglog>=1.4
  ${geometry_EXTRA_DEPENDS}
  DEFINITIONS ${geometry_EXTRA_DEFINITIONS})

set(geometry_SOURCES
  pointcloud.cpp
  volume.cpp
  parse-obj.cpp

  polygon.hpp
  polyclip.cpp

  triclip.hpp

  mesh.hpp mesh.cpp
  meshop.hpp meshop.cpp

  neighbors.hpp

  smoothcloud.hpp smoothcloud.cpp

  mesh-voxelizer.hpp mesh-voxelizer.cpp

  volumeop.hpp volumeop.cpp
  detail/volume.mcubes.hpp detail/volume.mcubes.cpp

  binmesh.hpp binmesh.cpp

  forsyth.hpp forsyth.cpp

  ${geometry_OPENMESH_SOURCES}
  ${geometry_CGAL_SOURCES}
  ${geometry_OPENCV_SOURCES}
)

add_library(geometry STATIC ${geometry_SOURCES})
buildsys_library(geometry)
target_link_libraries(geometry ${MODULE_LIBRARIES})
buildsys_target_compile_definitions(geometry ${MODULE_DEFINITIONS})

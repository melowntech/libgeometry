# bump version here
set(geometry_VERSION 1.11)

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

    volume.cpp
    volumeop.hpp volumeop.cpp
    )

  list(APPEND geometry_EXTRA_DEPENDS OpenCV)

  # extra sources that need imgproc
  if(MODULE_imgproc_FOUND)
    message(STATUS "geometry: compiling in imgproc support")
    list(APPEND geometry_EXTRA_DEFINITIONS GEOMETRY_HAS_IMGPROC=1)

    list(APPEND geometry_OPENCV_SOURCES
      mesh-voxelizer.hpp mesh-voxelizer.cpp
      )

    list(APPEND geometry_EXTRA_DEPENDS imgproc=>1.23)
  else()
    message(STATUS "geometry: compiling without imgproc support")
  endif()
else()
  message(STATUS "geometry: compiling without OpenCV support")
endif()


define_module(LIBRARY geometry=${geometry_VERSION}
  DEPENDS math>=1.4 utility>=1.11 dbglog>=1.4 gts glib
  ${geometry_EXTRA_DEPENDS}
  DEFINITIONS ${geometry_EXTRA_DEFINITIONS})

set(geometry_SOURCES
  pointcloud.cpp
  parse-obj.cpp

  polygon.hpp
  polyclip.cpp

  triclip.hpp

  mesh.hpp mesh.cpp
  meshop.hpp meshop.cpp

  neighbors.hpp

  smoothcloud.hpp smoothcloud.cpp


  detail/volume.mcubes.hpp detail/volume.mcubes.cpp

  binmesh.hpp binmesh.cpp

  forsyth.hpp forsyth.cpp

  nonconvexclip.hpp nonconvexclip.cpp
  triangulate.hpp triangulate.cpp
  meshop-gts.cpp

  ${geometry_OPENMESH_SOURCES}
  ${geometry_CGAL_SOURCES}
  ${geometry_OPENCV_SOURCES}
)

add_library(geometry STATIC ${geometry_SOURCES})
buildsys_library(geometry)
target_link_libraries(geometry ${MODULE_LIBRARIES})
buildsys_target_compile_definitions(geometry ${MODULE_DEFINITIONS})

# compatibility tweaking
if (${Boost_VERSION} LESS 105200)
  # Boost.Polygon before version 1.52 lacked Voronoi diagram support
  target_include_directories(geometry PRIVATE compat)
endif()

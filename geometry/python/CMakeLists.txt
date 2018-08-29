define_module(LIBRARY pygeometry=${MODULE_geo_VERSION}
  DEPENDS
  geometry=${MODULE_geometry_VERSION}
  Boost_THREAD Boost_SYSTEM Boost_PYTHON PYTHONLIBS)

set(pygeometry_SOURCES
  geometrymodule.cpp
  )

add_library(pygeometry STATIC ${pygeometry_SOURCES})
buildsys_library(pygeometry)
target_link_libraries(pygeometry ${MODULE_LIBRARIES})

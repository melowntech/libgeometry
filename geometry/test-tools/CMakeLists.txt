# test program
define_module(BINARY geometry-load-obj
  DEPENDS
  geometry
)

set(geometry-load-obj_SOURCES
  obj.cpp
  )

add_executable(geometry-load-obj ${geometry-load-obj_SOURCES})
target_link_libraries(geometry-load-obj ${MODULE_LIBRARIES})
target_compile_definitions(geometry-load-obj PRIVATE ${MODULE_DEFINITIONS})
buildsys_binary(geometry-load-obj)

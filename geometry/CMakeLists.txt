set(geometry_SOURCES
	pointcloud.cpp
	volume.cpp
    parse-obj.cpp
)

add_library(geometry STATIC ${geometry_SOURCES})

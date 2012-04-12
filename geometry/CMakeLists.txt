
set(geometry_SOURCES
	pointcloud.cpp
	volume.cpp
)

# This magic either creates (static) library or assigns all sources to variable
# passed in ${Add_Subdirectory_And_Grab_Sources_Var}.
#
# Used to fetch sources when static library is not an option (e.g. when building
# shared library containing object files compiled from these sources,
# (like python extension)).
if(Add_Subdirectory_And_Grab_Sources_Var)
  set(${Add_Subdirectory_And_Grab_Sources_Var}
    ${${Add_Subdirectory_And_Grab_Sources_Var}};${math_SOURCES}
    PARENT_SCOPE)
else()
  add_library(geometry STATIC ${geometry_SOURCES})
endif()

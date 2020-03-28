cmake_minimum_required (VERSION 3.1.3)
project (viewer_mpas_convertor)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (NetCDF REQUIRED)
include_directories (${NETCDF_INCLUDE_DIR})

include_directories ( include src/io src/utils)

# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable (main src/main.cpp src/io/mpaso.cpp src/io/particles.cpp src/io/segments.cpp src/utils/misc.cpp)
target_link_libraries (main ${NETCDF_LIBRARIES})



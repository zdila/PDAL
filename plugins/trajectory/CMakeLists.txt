project (trajectory-reconstruction)

cmake_minimum_required (VERSION 3.1.0) # This version was released 2014-12-15
set (PACKAGE_MAJOR_VERSION 1)
set (PACKAGE_MINOR_VERSION 0)
set (PACKAGE_PATCH_VERSION 0)

set (${PROJECT_NAME}_MAJOR_VERSION ${PACKAGE_MAJOR_VERSION})
set (${PROJECT_NAME}_MINOR_VERSION ${PACKAGE_MINOR_VERSION})
set (${PROJECT_NAME}_PATCH_VERSION ${PACKAGE_PATCH_VERSION})
set (PACKAGE_VERSION
  ${PACKAGE_MAJOR_VERSION}.${PACKAGE_MINOR_VERSION}.${PACKAGE_PATCH_VERSION})
set (${PROJECT_NAME}_VERSION ${PACKAGE_VERSION})
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  # Set a default build type for single-configuration cmake generators
  # if no build type is set.
  set (CMAKE_BUILD_TYPE Release)
endif ()
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (Eigen3 REQUIRED)
set (Eigen3_LIBRARIES Eigen3::Eigen)

find_package (Ceres REQUIRED)
find_package (GeographicLib REQUIRED)
find_package (PDAL 2.0 REQUIRED CONFIG)

if (MSVC)
  string (REGEX REPLACE "/W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
else ()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif ()

set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set (CMAKE_INSTALL_RPATH "\$ORIGIN/../lib${LIB_SUFFIX}")

set (PROJECT_LIBRARY_LIST LidarTrajectory)

add_subdirectory (src)
add_subdirectory (test)
add_subdirectory (filters)
add_subdirectory (resources)
add_subdirectory (cmake)
add_subdirectory (sandbox)

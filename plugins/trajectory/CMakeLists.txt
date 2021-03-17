find_package (Eigen3 REQUIRED)
set (Eigen3_LIBRARIES Eigen3::Eigen)
find_package (Ceres REQUIRED)
find_package (GeographicLib REQUIRED)

set (PROJECT_LIBRARY_LIST LidarTrajectory)

add_subdirectory (src)
add_subdirectory (filters)
add_subdirectory (resources)
add_subdirectory (cmake)
#add_subdirectory (sandbox)

if (WITH_TESTS)
    add_subdirectory (test)
endif()

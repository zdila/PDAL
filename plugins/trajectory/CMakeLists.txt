find_package (Eigen3 REQUIRED)
set (Eigen3_LIBRARIES Eigen3::Eigen)
find_package (Ceres REQUIRED)
find_package (GeographicLib REQUIRED)

set (PROJECT_LIBRARY_LIST LidarTrajectory)

add_subdirectory (src)
#add_subdirectory (filters)

PDAL_ADD_PLUGIN(filter_libname filter sritrajectory
    FILES
        filters/SRITrajectory.cpp
    LINK_WITH
        ${PROJECT_LIBRARY_LIST}
	${PDAL_LIBRARIES}
    )

add_subdirectory (resources)
add_subdirectory (cmake)
#add_subdirectory (sandbox)

# Currently not a standard test
#if (WITH_TESTS)
#    #add_subdirectory (test)
#    PDAL_ADD_TEST(pdal_filters_trajectory_test
#        FILES
#	    test/testtraj.cpp
#	LINK_WITH
#            ${PROJECT_LIBRARY_LIST}
#	)
#endif()

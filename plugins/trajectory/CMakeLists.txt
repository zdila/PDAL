find_package (Eigen3)
find_package (Ceres)
find_package (GeographicLib)

if ((NOT Ceres_FOUND) OR (NOT GeographicLib_FOUND) OR (NOT Eigen3_FOUND))
	message(WARNING "Ceres, GeographicLib or Eigen3 not found. Skipping trajectory filter plugin.")
else()
    set (Eigen3_LIBRARIES Eigen3::Eigen)
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
endif()

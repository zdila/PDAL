#
# Trajectory plugin CMake configuration
#
include(${PDAL_CMAKE_DIR}/ceres.cmake)
if (NOT CERES_FOUND)
    message(FATAL_ERROR "Can't find Ceres support required to build.")
endif()

#
# ScanAngleTrajectory Filter
#
PDAL_ADD_PLUGIN(filter_libname filter scanangletrajectory
    FILES
        filters/ScanAngleTrajectoryFilter.cpp
)

#
# ReturnInfoTrajectory Filter
#
PDAL_ADD_PLUGIN(filter_libname filter returninfotrajectory
    FILES
        filters/ReturnInfoTrajectoryFilter.cpp
    LINK_WITH
        ${CERES_LIBRARIES}
    INCLUDES
        ${CERES_INCLUDE_DIRS}
)

set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set (CMAKE_INSTALL_RPATH "\$ORIGIN/../lib${LIB_SUFFIX}")

add_library(pdal_plugin_filter_sritrajectory SHARED
  SRITrajectory.cpp SRITrajectory.hpp)

target_link_libraries(pdal_plugin_filter_sritrajectory
  PRIVATE ${PROJECT_LIBRARY_LIST} ${PDAL_LIBRARIES})

# install lib
install(FILES SRITrajectory.hpp DESTINATION include)
install(TARGETS
  pdal_plugin_filter_sritrajectory
  EXPORT targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})

# install (TARGETS qtm-tester EXPORT targets DESTINATION bin)

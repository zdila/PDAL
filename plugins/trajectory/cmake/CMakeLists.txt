#############################################################################
# Config file generation and installation
#############################################################################

# Set where cmake will install the proj-config.cmake file
# ${CMAKE_INSTALL_PREFIX}/${INSTALL_CMAKE_DIR} and @VISIONPROJ_ROOT_DIR@
# is the relative path to the root from there.
if (NOT WIN32)
  set(INSTALL_CMAKE_DIR "share/cmake/${PROJECT_NAME}")
  set (PROJECT_ROOT_DIR "../../..")
else ()
  set(INSTALL_CMAKE_DIR "cmake")
  set (PROJECT_ROOT_DIR "..")
endif ()

# Now create the project-config files using the .in templates in
# VisionProj/cmake
# @only means only expand @foo@ from the config template
configure_file (project-config.cmake.in
  "${PROJECT_NAME_LOWER}-config.cmake" @ONLY)
configure_file (project-config-version.cmake.in
  "${PROJECT_NAME_LOWER}-config-version.cmake" @ONLY)

# Install the project-config files so that other modules can find this project
install (FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")

# Make information about the cmake targets (the library)
# available.
#install (EXPORT targets
#  FILE "${PROJECT_NAME_LOWER}-targets.cmake"
#  DESTINATION "${INSTALL_CMAKE_DIR}")

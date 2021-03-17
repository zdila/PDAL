add_executable (traj-ceres traj-ceres.cpp lasertraj.hpp lasertraj.cpp)
target_include_directories (traj-ceres PRIVATE ${CERES_INCLUDE_DIRS})
target_link_libraries (traj-ceres
  PRIVATE ${Eigen3_LIBRARIES} ${CERES_LIBRARIES})

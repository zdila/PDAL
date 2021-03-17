import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

#actual_traj_file = 'D:/Dataset/preston_trajectory_data/sri_data/uh_campus/sbet_047_IGS08-UTM15N-Ellipsoid.txt'
actual_traj_file = 'D:/Dataset/preston_trajectory_data/sri_data/ak_sitka/apps_final_ATLANS-20160503_NAD83-UTM8N-Geoid12B.txt'

# pdal-sritrajectory result
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-test-scripts\\Results\\C2_L2_attitudefull_pitch2.txt'
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-test-scripts\\Results\\ak_sitka_1_attitude.txt'
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-weight-test\\Results\\C2_L2_attitude_multiweight_0001_scanweight_1.txt'
est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-weight-test\\Results\\ak_sitka_1_attitude_multiweight_1_scanweight_0001.txt'
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-weight-test\\Results\\ak_sitka_1_attitude_multiweight_0001_scanweight_1.txt'

# python result
#est_traj_file = 'D:/Dataset/preston_trajectory_data/sri_data/ak_sitka/ak_sitka_1_sorted_EstimatedTrajectory.txt'

# pdal-scanangle2 result
#est_traj_file = 'D:/Dataset/preston_trajectory_data/pdal-scanangle2-test-scripts/Results/ak_sitka_4_estimated_traj.txt'
#est_traj_file = 'D:/Dataset/preston_trajectory_data/pdal-scanangle2-test-scripts/Results/C2_L2_estimated_traj.txt'

# result from preston
#est_traj_file = 'D:/Dataset/preston_trajectory_data/pdal-scanangle2-test-scripts/Results/ak_sitka_1_estimated_traj.txt'


# actual_traj_file = 'D:/OneDrive - University Of Houston/2018-10-25 - ERDC BAA/ActualWork/sri/sri_data/ak_sitka/apps_final_ATLANS-20160503_NAD83-UTM8N-Geoid12B.txt'
# est_traj_file = 'D:/OneDrive - University Of Houston/2018-10-25 - ERDC BAA/ActualWork/sri/sri_data/ak_sitka/autoclass - Scanner 1 - 160503_011744_VQ480i - originalpoints_4-method2-EstimatedTrajectory.txt'

# actual_traj_file = 'F:/UH/sbet_047_IGS08-UTM15N-Ellipsoid.txt'
# est_traj_file = 'F:/UH/selective_ta5_mda15_j8.txt'

actual_traj = np.loadtxt(actual_traj_file, delimiter=',', skiprows=1)
est_traj = np.loadtxt(est_traj_file, delimiter=' ', skiprows=0)

# Trim actual trajectory to time bounds of estimated trajectory
min_t = np.min(est_traj[:,0])
max_t = np.max(est_traj[:,0])
mask = np.logical_and(actual_traj[:,0] >= min_t, actual_traj[:,0] <= max_t)
actual_traj = actual_traj[mask,:]

# "GpsTime","Y","X","Z","Roll","Pitch","Azimuth"
# actual_traj_x = actual_traj[:,2]
# actual_traj_y = actual_traj[:,1]
# actual_traj_z = actual_traj[:,3]

# Interpolate actual trajectory at estimated trajectory times for comparison
actual_interp_x = np.interp(est_traj[:,0], actual_traj[:,0], actual_traj[:,2])
actual_interp_y = np.interp(est_traj[:,0], actual_traj[:,0], actual_traj[:,1])
actual_interp_z = np.interp(est_traj[:,0], actual_traj[:,0], actual_traj[:,3])

min_x = min(actual_interp_x)
max_x = max(actual_interp_x)
min_y = min(actual_interp_y)
max_y = max(actual_interp_y)
min_z = min(actual_interp_z)
max_z = max(actual_interp_z)


# actual and est
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(actual_interp_x, actual_interp_y, actual_interp_z, 'gray', label='actual')
ax.plot3D(est_traj[:,1], est_traj[:,2], est_traj[:,3], 'red', label='est')
ax.set_xlim3d(min_x, max_x)
ax.set_ylim3d(min_y, max_y)
ax.set_zlim3d(min_z, max_z)
ax.legend()
plt.show()

# # actual
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot3D(actual_interp_x, actual_interp_y, actual_interp_z, 'gray')
# #ax.scatter3D(actual_interp_x, actual_interp_y, actual_interp_z, 'Greens')
# #ax.scatter3D(actual_traj[:,2], actual_traj[:,1], actual_traj[:,3], cmap='Greens');
# ax.set_xlim3d(min_x, max_x)
# ax.set_ylim3d(min_y, max_y)
# ax.set_zlim3d(min_z, max_z)
# plt.show()
#
# # est
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# # Data for three-dimensional scattered points
# #ax.scatter3D(est_traj[:,1], est_traj[:,2], est_traj[:,3], cmap='Greens');
# ax.plot3D(est_traj[:,1], est_traj[:,2], est_traj[:,3], 'red')
# ax.set_xlim3d(min_x, max_x)
# ax.set_ylim3d(min_y, max_y)
# ax.set_zlim3d(min_z, max_z)
# plt.show()

dummy=0
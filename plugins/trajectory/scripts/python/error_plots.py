import numpy as np
import matplotlib.pyplot as plt

use_360 = False

#actual_traj_file = 'D:/Dataset/preston_trajectory_data/sri_data/uh_campus/sbet_047_IGS08-UTM15N-Ellipsoid.txt'
actual_traj_file = 'D:/Dataset/preston_trajectory_data/sri_data/ak_sitka/apps_final_ATLANS-20160503_NAD83-UTM8N-Geoid12B.txt'

# pdal-sritrajectory result
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-weight-test\\Results\\C2_L2_attitude_multiweight_1_scanweight_0001.txt'
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-weight-test\\Results\\C2_L2_attitude_multiweight_0001_scanweight_1.txt'
est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-weight-test\\Results\\ak_sitka_1_attitude_multiweight_1_scanweight_0001_dtr_0005_dts_0005_pitch_nan.txt'
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-weight-test\\Results\\ak_sitka_1_attitude_multiweight_0001_scanweight_1.txt'

#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-test-scripts\\Results\\C2_L2_attitude.txt'
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-test-scripts\\Results\\C2_L2_attitudefull_pitch0.txt'
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-test-scripts\\Results\\ak_sitka_1_attitudefull_pitch0.txt'
#est_traj_file = 'D:\\Dataset\\preston_trajectory_data\\pdal-sritrajectory-test-scripts\\Results\\ak_sitka_1_attitude.txt'

# python result
#est_traj_file = 'D:/Dataset/preston_trajectory_data/sri_data/ak_sitka/ak_sitka_1_sorted_EstimatedTrajectory.txt'

# pdal-scanangle2 result
#est_traj_file = 'D:/Dataset/preston_trajectory_data/pdal-scanangle2-test-scripts/Results/ak_sitka_1_estimated_traj.txt'
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
# Interpolate actual trajectory at estimated trajectory times for comparison
actual_interp_x = np.interp(est_traj[:,0], actual_traj[:,0], actual_traj[:,2])
actual_interp_y = np.interp(est_traj[:,0], actual_traj[:,0], actual_traj[:,1])
actual_interp_z = np.interp(est_traj[:,0], actual_traj[:,0], actual_traj[:,3])

actual_interp_roll = np.interp(est_traj[:,0], actual_traj[:,0], actual_traj[:,4])
actual_interp_pitch = np.interp(est_traj[:,0], actual_traj[:,0], actual_traj[:,5])

# **this should be corrected: how to do rotation interpolation when you don't know which direction it's rotating?
actual_interp_yaw = np.interp(est_traj[:,0], actual_traj[:,0], actual_traj[:,6])

# errors
num_data_points = len(est_traj)
print("num_data_points: {}".format(num_data_points))

x_error = abs(actual_interp_x-est_traj[:,1])
x_error_max = max(x_error)
x_error_mean = np.mean(x_error)
print("x_error_max: {}".format(x_error_max))
print("x_error_avg: {}".format(x_error_mean))

y_error = abs(actual_interp_y-est_traj[:,2])
y_error_max = max(y_error)
y_error_mean = np.mean(y_error)
print("y_error_max: {}".format(y_error_max))
print("y_error_avg: {}".format(y_error_mean))

z_error = abs(actual_interp_z-est_traj[:,3])
z_error_max = max(z_error)
z_error_mean = np.mean(z_error)
print("z_error_max: {}".format(z_error_max))
print("z_error_avg: {}".format(z_error_mean))

if (est_traj.shape[1]>=6 ):
    # angles
    est_yaw = est_traj[:, 6]
    if use_360:
        est_yaw = (est_yaw % 360 + 360) % 360

    yaw_error = actual_interp_yaw-est_yaw
    yaw_error = (yaw_error + 180) % 360 - 180
    yaw_error_abs = abs(yaw_error)
    yaw_error_max = max(yaw_error_abs)
    yaw_error_mean = np.mean(yaw_error_abs)
    print("yaw_error_max: {}".format(yaw_error_max))
    print("yaw_error_avg: {}".format(yaw_error_mean))
else:
    print("yaw_error_max: n/a")
    print("yaw_error_avg: n/a")

if (est_traj.shape[1]>=6 ):
    est_pitch = est_traj[:,5]
    if use_360:
        est_pitch = (est_pitch % 360 + 360) % 360
    pitch_error = actual_interp_pitch-est_pitch
    pitch_error = (pitch_error + 180) % 360 - 180
    pitch_error_abs = abs(pitch_error)
    pitch_error_max = max(pitch_error_abs)
    pitch_error_mean = np.mean(pitch_error_abs)
    print("pitch_error_max: {}".format(pitch_error_max))
    print("pitch_error_avg: {}".format(pitch_error_mean))
else:
    print("pitch_error_max: n/a")
    print("pitch_error_avg: n/a")

traj_diff = np.diff(est_traj, axis=0);
dxdt = np.divide(traj_diff[:,1], traj_diff[:,0])
dydt = np.divide(traj_diff[:,2], traj_diff[:,0])
#headings = np.mod(-np.arctan2(traj_diff[:,1], traj_diff[:,0])*180/np.pi+180, 360)    # 0 to 360
headings = np.arctan2(dxdt, dydt)*180/np.pi                   # -180 to 180
headings = np.append(headings[0], headings)

if use_360:
    headings =  (headings + 360) % 360

heading_error = actual_interp_yaw-headings[:]
heading_error = (heading_error + 180) % 360 - 180   # smaller angle
heading_error_abs = abs(heading_error)

heading_error_max = max(heading_error_abs)
heading_error_mean = np.mean(heading_error_abs)
print("heading_error_max: {}".format(heading_error_max))
print("heading_error_avg: {}".format(heading_error_mean))

# Trajectory heights, height differences, roll/pitch/heading from actual trajectory
fig, axs =plt.subplots(4)
plt.subplots_adjust(hspace=0.5)
axs[0].plot(actual_traj[:,0], actual_traj[:,3], est_traj[:,0], est_traj[:,3])
axs[0].legend(('Actual Trajectory','Estimated Trajectory'))
axs[0].set(xlabel='Time (s)', ylabel='Height (m)', title='Trajectory Height')
axs[0].grid(True)
axs[1].plot(est_traj[:,0], actual_interp_z-est_traj[:,3])
axs[1].set(xlabel='Time (s)', ylabel='Height Difference(m)',title='Height Difference (Actual - Estimated)')
axs[1].grid(True)
axs[2].plot(actual_traj[:,0], actual_traj[:,4], 'g-', actual_traj[:,0], actual_traj[:,5], 'c-')
axs[2].legend(('Roll', 'Pitch'))
axs[2].set(xlabel='Time (s)', ylabel='Degrees', title='Actual Trajectory Roll and Pitch')
axs[2].grid(True)
axs[3].plot(actual_traj[:,0], actual_traj[:,6], 'g-')
axs[3].legend(('Actual Trajectory Heading'))
axs[3].set(xlabel='Time (s)', ylabel='Degrees', title='Actual Trajectory Heading')
axs[3].grid(True)
figure = plt.gcf()
figure.set_size_inches(18,9)
plt.show()
# plt.savefig("{}_{}_TrajectoryHeightDifference_{}_seconds_every_{} seconds.png".format(Area,flight,sample,every),bbox_inches='tight')

# Trajectory X and Y differences, roll/pitch/heading from actual trajectory
fig, axs =plt.subplots(4)
plt.subplots_adjust(hspace=0.5)
axs[0].plot(est_traj[:,0], actual_interp_x-est_traj[:,1])
axs[0].legend(('X Difference'))
axs[0].set(xlabel='Time (s)', ylabel='X Difference (m)', title='X Difference (Actual - Estimated)')
axs[0].grid(True)
axs[1].plot(est_traj[:,0], actual_interp_y-est_traj[:,2])
axs[1].legend(('Y Difference'))
axs[1].set(xlabel='Time (s)', ylabel='Y Difference (m)', title='Y Difference (Actual - Estimated)')
axs[1].grid(True)
axs[2].plot(actual_traj[:,0], actual_traj[:,4], 'g-',actual_traj[:,0], actual_traj[:,5], 'c-')
axs[2].legend(('Roll','Pitch'))
axs[2].set(xlabel='Time (s)', ylabel='Degrees', title='Actual Trajectory Roll and Pitch')
axs[2].grid(True)
axs[3].plot(actual_traj[:,0], actual_traj[:,6], 'g-')
axs[3].legend(('Heading'))
axs[3].set(xlabel='Time (s)', ylabel='Degrees', title='Actual Trajectory Heading')
axs[3].grid(True)
figure = plt.gcf()
figure.set_size_inches(18,9)
plt.show()
# plt.savefig("{}_{}_TrajectoryXYDifference_{}_seconds_every_{}_seconds.png".format(Area,flight,sample,every),bbox_inches='tight')

# yaw(Azimuth)/pitch estimation and difference
fig, axs =plt.subplots(4)
plt.subplots_adjust(hspace=0.5)
if (est_traj.shape[1]>=6 ):
    axs[0].plot(est_traj[:,0], yaw_error)
    axs[0].legend(('Yaw Difference'))
    axs[0].set(xlabel='Time (s)', ylabel='Yaw Difference (deg)', title='Yaw Difference (Actual - Estimated)')
    axs[0].grid(True)
else:
    #axs[0].plot(est_traj[:,0], actual_interp_pitch-est_traj[:,5])
    axs[0].plot(est_traj[:,0], heading_error)
    axs[0].legend(('Heading from XYZ Difference'))
    axs[0].set(xlabel='Time (s)', ylabel='Heading from XYZ Difference (deg)', title='Heading from XYZ Difference (Actual - Computed)')
    axs[0].grid(True)

if (est_traj.shape[1]>=5 ):
    axs[1].plot(est_traj[:,0], pitch_error)
    axs[1].legend(('Pitch Difference'))
    axs[1].set(xlabel='Time (s)', ylabel='Pitch Difference (deg)', title='Pitch Difference (Actual - Estimated)')
    axs[1].grid(True)
axs[2].plot(actual_traj[:,0], actual_traj[:,4], 'g-',actual_traj[:,0], actual_traj[:,5], 'c-')
axs[2].legend(('Roll','Pitch'))
axs[2].set(xlabel='Time (s)', ylabel='Degrees', title='Actual Trajectory Roll and Pitch')
axs[2].grid(True)


axs[3].plot(actual_traj[:,0], actual_traj[:,6], 'g-', est_traj[:,0], headings[:], 'c-')

if (est_traj.shape[1] >= 6):
    axs[3].plot(est_traj[:,0], est_yaw, 'r-')
    axs[3].legend(('Heading_gt', 'Heading_computed', 'Yaw estimated'))
else:
    axs[3].legend(('Heading_gt', 'Heading_computed'))

axs[3].set(xlabel='Time (s)', ylabel='Degrees', title='Actual Trajectory Heading')
axs[3].grid(True)
figure = plt.gcf()
figure.set_size_inches(18,9)
plt.show()
# plt.savefig("{}_{}_TrajectoryXYDifference_{}_seconds_every_{}_seconds.png".format(Area,flight,sample,every),bbox_inches='tight')



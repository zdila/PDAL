# [Trajectory estimation from multiple-return LIDAR pulses](README.md)

## [Abstract](01-abstract.md)

## Introduction

LIDAR data sets, typically provided in the form of `LAS` files, often
do not contain information about on the location of the sensor
platform as a function of time.  For datasets with include the GPS
time for each return, it is possible to identify the multiple returns
originating from a given LIDAR pulse.  By combining the data for
multiple pulses emitted in a short time, it is possible to
"triangulate" for the position of the sensor.  This idea was proposed
by Gatziolis and McGaughey (2019) who showed how to obtain a full
sensor trajectory.

Here we reformulate this problem with a view to obtaining a more
accurate trajectory.  The trajectory is modeled as a cubic spline fit.
Such a fit independently fits the $`x`$, $`y`$, and $`z`$ components
of $`\mathbf R(t)`$.  The *unknowns* in this model are the parameters
specifying the cubic splines.  The *knowns* are the positions (and
times) of the multiple returns from individual LIDAR pulses.  The
optimization problem is then to adjust parameters specifying the
trajectory to minimize the RMS error between the returns and a ray
drawn from the sensor position to the mean position of the return.
This is a rather complex nonlinear optimization problem.  Fortunately,
it is one that is easily handled by the the software library Ceres
Solver.

The errors in this problem are primarily quantization errors in the
positions of the returns.  In the normal post-processing of a LIDAR
collect, the positions of the returns are computed from IMU data from
the sensor platform, the scan angle of the LIDAR swee, and timing
information for the returns.  If these positions were recorded
accurately, it would be possible to determine a precise ray for a given
pulse (with 2 or more returns) and sensor position could be
accurately triangulated from the rays nearly simultaneous pulses.
However, the default precision of the return positions in a `LAS` file
is $`0.01\,\text{m}`$.  If the two returns are $`5\,\text{m}`$ apart
and the sensor is flying $`1000\,\text{m}`$ above the ground then,
then the possible rays consistent with the return data span
$`2\,\text{m}`$ at the altitude of the sensor.  The uncertainty in
the triangulation with such rays is increased because of ill
conditioned triangles (leading chiefly to a large uncertainty in the
height of the sensor).

## [Fixed sensor](03-fixed.md)

## [The trajectory computation](04-trad.md)

## [References](09-refs.md)

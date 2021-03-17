# [Trajectory estimation from multiple-return LIDAR pulses](README.md)

## [Abstract](01-abstract.md)

## [Introduction](02-intro.md)

## [Fixed sensor](03-fixed.md)

## The trajectory computation

The discussion above solves for the sensor position at a single instant
of time.  Of course, the sensor position is typically moving and it is
convenient to model the motion as a cubic spline.  One approach would be
to perform a series of fixed sensor calculations, e.g., at
$`0.01\,\text s`$ intervals including for each calculation 10 pulses
sampled at $`0.001\,\text s`$ intervals and then to fit a spline to
the resulting positions.

This approach has the drawback that some of the positions may be better
approximated than others and the spline fit should respect this.  This
could be achieved by assigning weights to the various position estimated
and this, essentially, is how Gatziolis and McGaughey addressed this
issue.  However this put another layer of complexity into the problem.

However, in the spirit of Ceres Solver, it make more sense to pose the
entire exercise as a single least-squares problem.  Let's start by
describing how to express a cubic spline.

### The cubic spline

A cubic spline is a piece-wise cubic polynomial function which in our
application we will use to approximate $`\mathbf R(t)`$.  Each
component of $`\mathbf R(t)`$ can be fit independently of the others.
So we only need to consider a cubic spline for a scalar function
$`f(t)`$ defined between $`T_0`$ and $`T_K = T_0 + K \,\Delta t`$.
The time interval divided in $`K`$ block of duration $`\Delta t`$,
with the $`k\text{th}`$ block consisting of the interval $`T_k \le t
< T_{k+1}`$ where $`T_k = T_0 + k \,\Delta t`$ and $`k \in [0,K)`$.
Between each block, $`T_k`$ for $`k \in (0, K)`$, we require that
$`f(T_k)`$, $`f'(T_k)`$, and $`f''(T_k)`$ are continuous.

We shall specify the cubic polynomial for the $`k\text{th}`$ block
by the values of $`f(t)`$ and $`g(t) = \Delta t\,f'(t)`$ at the block
boundaries.  It is convenient to introduced a scaled centered time
variable the block $`\tau = (t - T_k)/\Delta t - \frac12`$.  At the
block boundaries, we have $`f = f_k`$ and $`g = g_k`$ at $`\tau =
-\frac12`$ and $`f = f_{k+1}`$ and $`g = g_{k+1}`$ at $`\tau =
\frac12`$.  It is now a simple matter, e.g., by using the algebra
system, Maxima, to find the polynomial satisfying the boundary
conditions

```math
f(\tau) = a_0 + a_1 \tau + a_2 \tau^2 + a_3 \tau^3,
```

where

```math
\begin{aligned}
 a_0 &= \textstyle\frac18 (4 f_+ - g_-),\\
 a_1 &= \textstyle\frac14 (6 f_- - g_+),\\
 a_2 &= \textstyle\frac12 g_-,\\
 a_3 &= -2 f_- + g_+,\\
 f_+ &= f_{k+1} + f_k, \\
 f_- &= f_{k+1} - f_k, \\
 g_+ &= g_{k+1} + g_k, \\
 g_- &= g_{k+1} - g_k.
\end{aligned}
```

By specifying the cubic polynomial by its values and derivatives at
the block boundaries, we ensure continuity of $`f(t)`$ and $`f'(t)`$.
The jump in the second derivative is $`j_k/\Delta t^2`$ where

```math
 j_k = 6 (f_{k+1}-f_{k-1}) - 2 (g_{k+1}+g_{k-1}) - 8 g_k.
```

We will then add $`j_k \approx 0`$ to the optimization problem.

### The optimization problem

We are now ready to set up the optimization problem for the entire
trajectory.  The *knowns* are $`t_i`$, $`\mathbf r_i`$, $`\mathbf
p_i`$, and $`d_i`$ for $`n`$ pulses.  These are the same as for the
fixed sensor problem with the addition of the time $`t_i`$ for each
pulse.  The *unknowns* are the sensor positions, $`R(t)`$ and
velocities, $`R'(t)`$, at the block boundaries $`t = T_k`$ for $`k \in
[0,K]`$.

The residue block for Ceres Solver now includes the time $`t_i`$.
Each pulse is assigned to a particular time block and the constructor
converts this time to the scaled time $`\tau`$.  The corresponding
function object now takes the trajectory position and velocities at
the block boundaries as input.  It then uses the stored value of
$`\tau`$ to evaluate the corresponding cubic polynomials for the 3
components of the sensor position $`R(t_i)`$.  The calculation then
proceeds as in the fixed-sensor case and returns a two-component
vector for the residue.

There is now also a new type of residue block to enforce the
continuity of the acceleration at the internal block boundaries.  The
function object takes $`R(T_{k-1})`$, $`R'(T_{k-1})`$, $`R'(T_k)`$,
$`R(T_{k+1})`$, and $`R'(T_{k+1})`$, and returns the jump in
$`R''(T_k)`$ which is proportional to $`j_k`$.

Overall there are $`6(K+1)`$ unknowns (the positions and velocities at
the block boundaries).  The number of equations is $`2n`$ for the
pulse residues plus $`3(K-1)`$ for the acceleration jump constraints.

We have tested this on flights lasting about a minute with $`\Delta t
= 0.1\,\text s`$ and sampling one multi-return pulse every
$`0.001\,\text s`$ from the LIDAR data (we select the pulse with the
largest distance between its first and last returns).  Even though
this involves a system of tens of thousands of equations, Ceres Solver
handles it without difficulty in a few seconds of CPU time.

## [References](09-refs.md)

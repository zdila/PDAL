# [Trajectory estimation from multiple-return LIDAR pulses](README.md)

## [Abstract](01-abstract.md)

## [Introduction](02-intro.md)

## Fixed sensor

In order to introduce the concepts, let us start first by assuming
that the sensor is fixed and emits $`n`$ multi-return pulses, indexed
by $`i \in [1,n]`$.  We shall only consider the first and last returns
(ignoring any intermediate returns).  We denote positions of the
returns by

```math
 \mathbf r_i^\plusmn = \mathbf r_i \plusmn d_i \mathbf p_i,
```

where superscripts $`\plusmn`$ denote first and last returns, $`\mathbf
p_i`$ is the unit vector in the direction from the last to the first
return, and $`d_i`$ is half the distance between the returns.

### The inverse method

The goal now is to determine the position $`\mathbf R`$ consistent
with these returns.  One approach is to consider the $`n`$ rays

```math
 \mathbf r_i + s_i \mathbf p_i
```

where the distance along the ray is parameterized by $`s_i`$ and to
solve the $`3n`$ equations

```math
 \mathbf r_i + s_i \mathbf p_i - \mathbf R = \mathbf h_i \approx 0
```

for the $`3+n`$ unknowns $`\mathbf R`$ and $`s_i`$.  This is an
overdetermined system for 2 or more pulses and then we can use
standard linear algebra methods to find the best solution in the
least-squares sense.

This is the approach used by Gatziolis and McGaughey who consider just
pairs of pulses $`n = 2`$.  The problem is that the resulting solution
for $`\mathbf R`$ is typically not the optimal solution for the
trajectory problem because the system of equations does not involve
the return separation $`d_i`$ so pulses with closely separated returns
and treated equally to pulses with widely separated returns.  In
reality the latter returns should be weighted more heavily.

Gatziolis and McGaughey address this problem by selecting an *optimal*
pair of returns based on the return separation and the angle between
the pulses.  This is based on a weighting function which needs to be
separately estimated.

We use a simplified version of this linear least squares problem to
find an initial trajectory for our method described below.  We use the
$`z`$ component of the residue equations to eliminate $`s_i`$ from the
system.  The equations are then

```math
\begin{aligned}
\biggl( r_{i,x} - \frac{p_{i,x}}{p_{i,z}} r_{i,z} \biggr) -
\biggl( R_x - \frac{p_{i,x}}{p_{i,z}} R_z \biggr) = h_{i,x} \approx 0, \\
\biggl( r_{i,y} - \frac{p_{i,y}}{p_{i,z}} r_{i,z} \biggr) -
\biggl( R_y - \frac{p_{i,y}}{p_{i,z}} R_z \biggr) = h_{i,y} \approx 0.
\end{aligned}
```

We can write this as the explicit overdetermined linear system

```math
\mathsf A \cdot \mathbf R = \mathbf B
```

where $`\mathsf A`$ is the $`2n \times 3`$ matrix

```math
\mathsf A = \begin{pmatrix}
\vdots & \vdots & \vdots \\
1 & 0 & - p_{i,x}/p_{i,z} \\
0 & 1 & - p_{i,y}/p_{i,z} \\
\vdots & \vdots & \vdots
\end{pmatrix}
```

(two rows for each of the $`n`$ pulses), $`\mathbf B`$ is the $`2n`$
column vector

```math
\mathbf B = \begin{pmatrix}
\vdots \\
r_{i,x} - (p_{i,x}/p_{i,z}) r_{i,z} \\
r_{i,y} - (p_{i,y}/p_{i,z}) r_{i,z} \\
\vdots
\end{pmatrix},
```

and $`\mathbf R`$ is the unknown sensor position to solve for.

This reduces the problem to $`2n`$ equations for $`3`$ unknowns.  In
this formulation we determine the horizontal plane $`z = R_z`$ in
which the rays are most tightly clustered.  This is *not* the same
problem as before; however with typical aerial LIDAR collects the two
solutions for $`R`$ will be reasonably close.  The difference is
immaterial in our application since this solution for $`R`$ is only
used as an initial estimate.

### The direct method

We term this of estimating $`\mathbf R`$ the *inverse* method,
because the rays are traced back from the returns to the sensor.  An
alternative is to trace the rays from $`\mathbf R`$ to the midpoint of
returns, the *direct* method.  Thus each ray is given by

```math
 \mathbf r_i + s_i \mathbf q_i,
```

where $`\mathbf q_i`$ is the unit vector from $`\mathbf r_i`$ to
$`\mathbf R`$.  The rays now all intersect at $`\mathbf R`$ and the
optimization problem is to solve

```math
  d_i \mathbf p_i - s_i \mathbf q_i = \mathbf e_i \approx 0.
```

Again we have $`3n`$ equations with $`3 + n`$ unknowns.  However the
quantities that are being minimized, $`\mathbf e_i`$, is the distance
between the ray and the given positions of the returns.  This method
now naturally gives more weight to widely separated returns and the
solution will similarly be more heavily governed by rays forming
well-conditioned triangles.

Incidentally, the ray from $`\mathbf R`$ to $`\mathbf r_i`$ passes
equally close to the first and last returns, so it is only necessary
to the minimize the distance to the first returns.

This system of equations is no longer linear, so it cannot be solved
by linear algebra techniques.  However, it is ideally suited for the
Ceres Solver package.  This finds the least-squares solution for
nonlinear optimization problems.  It also features

* Automatic determination of the Jacobian needed to find the
  solutions.  This is achieved by writing the formulas in standard
  notation but with the variables having a C++ type `Jet` which
  combines a quantity and its derivative and, through overloaded
  operators and functions, follows all the standard rules of
  differentiation.

* A robust optimization.  A standard problem of least-squares methods
  is that outliers in the data can skew the solution away from the
  "right" one.  Ceres Solver includes a variety of "loss functions"
  which cause the effect of errors in the equations to fall off past
  some threshold.  For example in this case, the threshold for the
  loss functions might be set to $`0.01\,\text{m}`$

### Simplifying the direct method

We can simplify the problem by observing that $`\mathbf e_i`$ is
minimized with $`s_i \approx d_i`$ and that the resulting $`\mathbf
e_i`$ then spans a two-dimensional space perpendicular to $`\mathbf
p_i`$.  Thus we can approximate the error $`\mathbf e_i`$ by
projecting the $`\mathbf q_i`$ onto the plane perpendicular to
$`\mathbf p_i`$ at the first return.  The first step is to convert to
a primed coordinate system with $`\mathbf r_i`$ at the origin and with
the $`\mathbf z'`$ axis parallel to $`\mathbf p_i`$.  This is achieved
by the rotation matrix

```math
 \mathsf M_i = \begin{pmatrix}
\displaystyle
\frac{p_{i,x}^2 p_{i,z} + p_{i,y}^2}{p_{i,x}^2+p_{i,y}^2} &
\displaystyle
\frac{-(1 - p_{i,z}) p_{i,x} p_{i,y}}{p_{i,x}^2+p_{i,y}^2} &
-p_{i,x} \\
\displaystyle
\frac{-(1 - p_{i,z}) p_{i,x} p_{i,y}}{p_{i,x}^2+p_{i,y}^2} &
\displaystyle
\frac{p_{i,x}^2 + p_{i,y}^2 p_{i,z}}{p_{i,x}^2+p_{i,y}^2} &
-p_{i,y} \\
p_{i,x} & p_{i,y} & p_{i,z}
\end{pmatrix}.
```

This matrix rotates the coordinate system about the axis $`\mathbf z
\times \mathbf p_i`$.  Applying this translation and rotation to the
sensor position gives

```math
 \mathbf R'_i = \mathsf M_i \cdot (\mathbf R - \mathbf r_i)
```

Finally we project $`\mathbf R'_i`$ onto the plane $`z' = d`$ which gives

```math
\mathbf e_i = \frac{d_i}{R'_{i,z}}
\begin{pmatrix}
R'_{i,x}\\
R'_{i,y}
\end{pmatrix}.
```

Now the number of unknowns is just $`3`$, the coordinates of $`\mathbf
R`$, and the number of equations is $`2n`$, $`\mathbf e_i \approx 0`$
for each two-component vector $`\mathbf e_i`$

Solving this least-squares problem with Ceres Solver entails writing a
C++ class implementing a "residue block".  The constructor for the
class takes the *knowns* for a particular pulse, i.e., $`\mathbf
r_i`$, $`\mathbf p_i`$, and $`d_i`$, and implements a function object
which accepts the *unknowns* $`\mathbf R`$ as input and returns the
residue $`\mathbf e_i`$.  This entails merely expressing the equations
above as computer code.  The problem is specified by $`n`$ such
residue blocks and an initial guess for $`\mathbf R`$ (obtained, for
example, by the inverse linear least squares problem).  Ceres Solver
repeatedly invokes the function objects while adjusting $`\mathbf R`$
to minimize $`\sum e_i^2`$.  Because of the automatic differentiation
built into Ceres Solver, it can compute the Jacobian for the problem
which says how each component of $`\mathbf e_i`$ changes as each
component of $`\mathbf R`$ is varied.  This allows Ceres Solver to
vary $`\mathbf R`$ in an optimal way in its search for the
least-squares solution.

## [The trajectory computation](04-traj.md)

## [References](09-refs.md)

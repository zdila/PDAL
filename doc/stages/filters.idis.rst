.. _filters.idis:

filters.idis
===============================================================================

The **Inverse Density Imporance Sampling Filter** subsamples the input
``PointView`` by first computing the inverse density of a k-neighborhood of
points where the density is defined as the summation of Euclidean distances to
each of the k-neighbors. The specified number of points (``count``) with the
lowest inverse density value are then appended to the output ``PointView``. See
[Groh2018]_ for further details.

.. embed::

Options
-------------------------------------------------------------------------------

count
  Desired number of output samples. [Default: 1000]

knn
  Number of neighbors used to compute density. [Default: 16]

.. include:: filter_opts.rst


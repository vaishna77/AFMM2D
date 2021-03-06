Benchmarks
^^^^^^^^^^

All the following benchmarks have been carried out on a 2.3 GHz Intel Core i5 processor with 8GB RAM (with OpenMP enabled, this is 4 threads), with g++-9 compiler and Eigen version 3.3.7. The compiler flags that were utilized are the same are those mentioned in the CMakeLists.txt file.

Time Taken vs Tolerance
~~~~~~~~~~~~~~~~~~~~~~~

These benchmarks were performed for size of the matrix :math:`N = 1000000`, with the size of the leaf node set to :math:`M = 100`.

+----------------+------------+---------+
|Tolerance       | Assembly(s)|MatVec(s)|
+================+============+=========+
|:math:`10^{-2}` | 19.8358    | 9.65095 |
+----------------+------------+---------+
|:math:`10^{-4}` | 49.7557    | 8.55117 |
+----------------+------------+---------+
|:math:`10^{-6}` | 70.4369    | 8.57725 |
+----------------+------------+---------+
|:math:`10^{-8}` | 126.125    | 12.9264 |
+----------------+------------+---------+
|:math:`10^{-10}`| 214.757    | 52.5742 |
+----------------+------------+---------+
|:math:`10^{-12}`| 419.34     | 47.0097 |
+----------------+------------+---------+
|:math:`10^{-14}`| 749.672    | 79.0786 |
+----------------+------------+---------+

.. image:: images/timeVsEpsilon.png
   :width: 600


Time Taken vs Size of Matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For these benchmarks, the leaf size was fixed at :math:`M = 100`, with tolerance set to :math:`10^{-12}`

+-----------------------+------------+------------+
|:math:`N`              | Assembly(s)|MatVec(s)   |
+=======================+============+============+
|:math:`5 \times 10^{3}`| 0.309096   | 0.040681   |
+-----------------------+------------+------------+
|:math:`10^{4}`         | 1.88374    | 0.062273   |
+-----------------------+------------+------------+
|:math:`5 \times 10^{4}`| 14.2385    | 0.47138    |
+-----------------------+------------+------------+
|:math:`10^{5}`         | 29.3395    | 1.20935    |
+-----------------------+------------+------------+
|:math:`5 \times 10^{5}`| 244.375    | 42.0531    |
+-----------------------+------------+------------+
|:math:`10^{6}`         | 419.34     | 47.0097    |
+-----------------------+------------+------------+


.. image:: images/timeVsN.png
   :width: 600

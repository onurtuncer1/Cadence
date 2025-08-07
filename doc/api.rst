------------------------------------------------------------------------------
.. Project: Cadence Math Library
.. Copyright (c) 2025, Onur Tuncer, PhD, Istanbul Technical University
..
.. SPDX-License-Identifier: MIT
.. License-Filename: LICENSE
------------------------------------------------------------------------------

.. _api:

Cadence C++ API Reference
=========================

This is the API reference for the Cadence Math Library. It includes templated classes for fixed-size vectors, matrices, and transformation utilities used in computational geometry, physics, robotics, and CAD.

----

.. doxygenclass:: Cadence::Vector
   :project: Cadence
   :members:
   :protected-members:
   :private-members:
   :undoc-members:
   :outline:
   :allow-dot-graphs:

.. doxygenclass:: Cadence::Transform
   :project: Cadence
   :members:
   :protected-members:
   :private-members:
   :undoc-members:
   :outline:
   :allow-dot-graphs:

.. doxygenclass:: Cadence::Matrix
   :project: Cadence
   :members:
   :protected-members:
   :private-members:
   :undoc-members:
   :outline:
   :allow-dot-graphs:

----

Helper Concepts and Type Aliases
================================

.. doxygentypedef:: Cadence::Vector2
   :project: Cadence

.. doxygentypedef:: Cadence::Vector3
   :project: Cadence

.. doxygentypedef:: Cadence::Vector4
   :project: Cadence

.. doxygentypedef:: Cadence::Matrix2x2
   :project: Cadence

.. doxygentypedef:: Cadence::Matrix3x3
   :project: Cadence

.. doxygentypedef:: Cadence::Matrix4x4
   :project: Cadence

.. doxygenconcept:: Cadence::MatrixExpression
   :project: Cadence

.. doxygenconcept:: Cadence::SymmetricPositiveDefinite
   :project: Cadence

----

Free Functions
==============

.. doxygenfunction:: Cadence::determinant(const Matrix<T, 1, 1>&)
   :project: Cadence

.. doxygenfunction:: Cadence::determinant(const Matrix<T, 2, 2>&)
   :project: Cadence

.. doxygenfunction:: Cadence::determinant(const Matrix<T, 3, 3>&)
   :project: Cadence

.. doxygenfunction:: Cadence::determinant(const Matrix<T, N, N>&)
   :project: Cadence

.. doxygenfunction:: Cadence::minor_matrix(const Matrix<T, N, N>&, std::size_t, std::size_t)
   :project: Cadence

.. doxygenfunction:: Cadence::is_symmetric(const Matrix<T, N, N>&, T)
   :project: Cadence

.. doxygenfunction:: Cadence::is_positive_definite(const Matrix<T, N, N>&, T)
   :project: Cadence

.. doxygenfunction:: Cadence::cholesky(const Matrix<T, N, N>&)
   :project: Cadence

----

.. note::

   This API is under active development. Contributions are welcome. For issues or feature requests, please visit the `Cadence GitHub repository <https://github.com/onurtuncer1/Cadence>`_.

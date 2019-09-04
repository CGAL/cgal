// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic2Triangulation2Concepts
\cgalConcept

The concept `Periodic_2Offset_2` describes a two-/dimensional integer vector with some specialized access functions and operations.

\cgalHasModel `CGAL::Periodic_2_offset_2`

\sa `Periodic_2TriangulationTraits_2`
\sa `Periodic_2DelaunayTriangulationTraits_2`

*/

class Periodic_2Offset_2
{
public:

/// \name Creation
/// @{

  /*!
  Default constructor.
  */
  Periodic_2Offset_2();

  /*!
  Constructs the offset (x,y).
  */
  Periodic_2Offset_2(int x, int y);

/// @}

/// \name Operations
/// @{

  /*!
  Return the vector sum of `this` and `o`.
  */
  Periodic_2Offset_2 operator+(const Periodic_2Offset_2 & o)
  const;

  /*!
  Return the vector difference of `this` and `o`.
  */
  Periodic_2Offset_2 operator-(const Periodic_2Offset_2 & o)
  const;

  /*!
  Return the negative vector of `this`.
  */
  Periodic_2Offset_2 operator-() const;

  /*!
  Add `o` to `this` using vector addition.
  */
  void operator+=(const Periodic_2Offset_2 & o)
  const;

  /*!
  Subtract `o` from `this` using vector subtraction.
  */
  void operator-=(const Periodic_2Offset_2 & o)
  const;

  /*!
  Return `true` if `o` and `this` represent the same vector.
  */
  bool operator==(const Periodic_2Offset_2 & o) const;

  /*!
  Return `true` if `o` and `this` do not represent the same
  vector.
  */
  bool operator!=(const Periodic_2Offset_2 & o) const;

  /*!
  Compare `this` and `o` lexicographically.
  */
  bool operator<(const Periodic_2Offset_2 & o) const;

/// @}

/// \name Access Functions
/// @{

  /*!
  Return the \f$ i\f$-th entry of `this`.
  \pre \f$ i\in\{0,1\}\f$
  */
  int operator[](int i);

  /*!
  Return the \f$ x\f$-entry of `this`.
  */
  int x() const;

  /*!
  Return the \f$ y\f$-entry of `this`.
  */
  int y() const;

  /*!
  Returns `true` if `this` is equal to (0,0).
  */
  bool is_null() const;

  /*!
  Returns `true` if `this` is equal to (0,0).
  */
  bool is_zero() const;

/// @}

}; /* end Periodic_2Offset_2 */


/*!
Inputs an offset from `is`.
\relates Periodic_2Offset_2
*/
istream& operator>>(istream & is, Periodic_2Offset_2 & off);

/*!
Outputs an offset from `os`.
\relates Periodic_2Offset_2
*/
ostream& operator<<(ostream & os, Periodic_2Offset_2 & off) const;

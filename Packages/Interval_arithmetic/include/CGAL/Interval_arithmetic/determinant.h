// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Interval_arithmetic/determinant.h
// package       : Interval_arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

// This is a file containing overloading functions of
// include/CGAL/determinant.h for the Interval_nt_advanced type.
// The most important one is sign_of_det_2x2() for the moment.
//
// Sylvain Pion.

template <>
inline
Sign
sign_of_determinant2x2( const Interval_nt_advanced& px,
			const Interval_nt_advanced& py,
                        const Interval_nt_advanced& qx,
			const Interval_nt_advanced& qy)
{
  // We have to determine the extreme point (the corners of the intervals that
  // define their cones wrt the origin).

  double p1x, p2x, p1y, p2y, q1x, q2x, q1y, q2y;
  unsigned int pc = 0, qc = 0;

// Sets p1x and p2x.

  if (py.inf() >= 0) // P > (Ox) // These tests must be strict for the filter.
  {
    p1x = px.sup();
    p2x = px.inf();
  }
  else if (py.sup() <= 0) // P < (Ox)
  {
    p1x = px.inf();
    p2x = px.sup();
    pc = 2;
  }
  else // P cuts (Ox)
  { // We decide to use the generic function here ?
    if (px.inf() >= 0) // P > (Oy)
      p1x = p2x = px.inf();
    else if (px.sup() <= 0) // P < (Oy)
      p1x = p2x = px.sup();
    else return // P contains 0 => sign unknown => throw exception.
  }

// Sets p1y and p2y.

  if (px.inf() >= 0) // P > (Oy)
  {
    p1y = py.inf();
    p2y = py.sup();
    if (pc!=0) pc = 3;
  }
  else if (px.sup() <= 0) // P < (Oy)
  {
    p1y = py.sup();
    p2y = py.inf();
    if (pc!=2) pc = 1;
  }
  else // P cuts (Oy)
  {
    if (py.inf() >= 0) // P > (Ox)
      p1y = p2y = py.inf();
    else if (py.sup() <= 0) // P < (Ox)    // Useless test, previously tested
      p1y = p2y = py.sup();
    else return // P contains 0.  Cannot happen, already tested.
  }

// Admit the same is done for q{12}{xy}.  (we can also already conclude from
// there sooner in some cases)

// What we now have to compare is the sign_of_det2x2(p1,q2) AND
// sign_of_det2x2(p2,q1), right ?


// What the default template version does :
// return static_cast<Sign>(static_cast<int>(compare( px*qy, py*qx)));
}


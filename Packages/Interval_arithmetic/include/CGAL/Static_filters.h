// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Static_filters.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================
 
#ifndef CGAL_STATIC_FILTERS_H
#define CGAL_STATIC_FILTERS_H

#include <CGAL/basic.h>

// Workaround for buggy compilers.
#ifdef CGAL_CFG_MATCHING_BUG_2
#  define CGAL_IA_CT double
#  define CGAL_IA_ET CGAL::MP_Float
#  define CGAL_IA_PROTECTED true
#  define CGAL_IA_CACHE No_Filter_Cache
#endif

#include <CGAL/Static_filters/Orientation_2.h>
#include <CGAL/Static_filters/Orientation_3.h>
#include <CGAL/Static_filters/Side_of_oriented_circle_2.h>
#include <CGAL/Static_filters/Side_of_oriented_sphere_3.h>
#include <CGAL/Static_filters/Coplanar_orientation_3.h>
#include <CGAL/Static_filters/Coplanar_side_of_bounded_circle_3.h>

// This traits class gathers optimized predicates written by hand, using
// a few steps of filtering.  It should work if the initial traits has
// cartesian coordinates which fit exactly in doubles.
//
// To allow pure static filters, the constant bound is computed by calling
// register_object() over each object (currently only Point_3) that some
// predicate of the traits will have to deal with.
//
// Note that we may extract more information with that : for example decide
// if all coordinates are fixed points values of a maximum number of bits,
// which guarantees the initial subtractions are exact...

CGAL_BEGIN_NAMESPACE

// Utility function to check a posteriori that a subtraction was performed
// without rounding error.
inline bool diff_was_exact(double a, double b, double ab)
{
    return ab+b == a && a-ab == b;
}

template < class K_base >
class Static_filters : public K_base
{
public :

  typedef typename K_base::Point_2 Point_2;
  typedef typename K_base::Point_3 Point_3;

  Static_filters()
    : max2x(0), max2y(0),
      max3x(0), max3y(0), max3z(0) {}

  typedef SF_Orientation_2<Point_2>                 Orientation_2;
  typedef SF_Orientation_3<Point_3>                 Orientation_3;
  typedef SF_Side_of_oriented_circle_2<Point_2>     Side_of_oriented_circle_2;
  typedef SF_Side_of_oriented_sphere_3<Point_3>     Side_of_oriented_sphere_3;
  typedef SF_Coplanar_orientation_3<Point_3, Orientation_2>
                                                    Coplanar_orientation_3;
  typedef SF_Side_of_bounded_circle_3<Point_3>
                                            Coplanar_side_of_bounded_circle_3;

  const Orientation_2 &
  orientation_2_object() const
  { return _orientation_2; }

  const Orientation_3 &
  orientation_3_object() const
  { return _orientation_3; }

  const Side_of_oriented_circle_2 &
  side_of_oriented_circle_2_object() const
  { return _side_of_oriented_circle_2; }

  const Side_of_oriented_sphere_3 &
  side_of_oriented_sphere_3_object() const
  { return _side_of_oriented_sphere_3; }

  const Coplanar_orientation_3 &
  coplanar_orientation_3_object() const
  { return _coplanar_orientation_3; }

  const Coplanar_side_of_bounded_circle_3 &
  coplanar_side_of_bounded_circle_3_object() const
  { return _coplanar_side_of_bounded_circle_3; }

  // These should not be const, but unfortunately Triangulation_?.geom_traits()
  // only give a const& access (should this be changed ?).
  // In the mean time, I put the data members mutable.

  void register_object(const Point_2 &p) const
  {
      bool redo = false;
      double dx = fabs(CGAL::to_double(p.x()));
      if (dx > max2x)
	  max2x = dx, redo = true;
      double dy = fabs(CGAL::to_double(p.y()));
      if (dy > max2y)
	  max2y = dy, redo = true;
      if (redo) {
          _orientation_2.update(max2x, max2y);
          _side_of_oriented_circle_2.update(max2x, max2y);
      }
  }

  void register_object(const Point_3 &p) const
  {
      bool redo = false;
      double dx = fabs(CGAL::to_double(p.x()));
      if (dx > max3x)
	  max3x = dx, redo = true;
      double dy = fabs(CGAL::to_double(p.y()));
      if (dy > max3y)
	  max3y = dy, redo = true;
      double dz = fabs(CGAL::to_double(p.z()));
      if (dx > max3z)
	  max3z = dz, redo = true;
      if (redo) {
          _orientation_3.update(max3x, max3y, max3z);
          _side_of_oriented_sphere_3.update(max3x, max3y, max3z);
	  _coplanar_orientation_3.update(max3x, max3y, max3z);
	  _coplanar_side_of_bounded_circle_3.update(max3x, max3y, max3z);
      }
  }

private:
  // Bounds on fabs() of the coordinates of the Point_2s and Point_3s.
  mutable double max2x, max2y;
  mutable double max3x, max3y, max3z;

  // A data member for each predicate.
  // Their state is related to the state of *this.
  // TODO : maybe they should be separate classes, with a pointer to the
  // kernel ?  That way we could copy them.
  mutable Orientation_2                     _orientation_2;
  mutable Orientation_3                     _orientation_3;
  mutable Side_of_oriented_circle_2         _side_of_oriented_circle_2;
  mutable Side_of_oriented_sphere_3         _side_of_oriented_sphere_3;
  mutable Coplanar_orientation_3            _coplanar_orientation_3;
  mutable Coplanar_side_of_bounded_circle_3 _coplanar_side_of_bounded_circle_3;
};

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_H

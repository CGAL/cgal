// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/predicates_on_points_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_PREDICATES_ON_POINTS_3_H
#define CGAL_CARTESIAN_PREDICATES_ON_POINTS_3_H

#include <CGAL/predicates/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
bool
x_equal(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().equal_x_3_object()(p, q);
}

template < class K >
inline
bool
y_equal(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().equal_y_3_object()(p, q);
}

template < class K >
inline
bool
z_equal(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().equal_z_3_object()(p, q);
}

template < class K >
inline
bool
less_x(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().less_x_3_object()(p, q);
}

template < class K >
inline
bool
less_y(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().less_y_3_object()(p, q);
}

template < class K >
inline
bool
less_z(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().less_z_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_x(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().compare_x_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_y(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().compare_y_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_z(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().compare_z_3_object()(p, q);
}

template < class K >
inline
bool
equal_xy(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().equal_xy_3_object()(p, q);
}

template < class K >
inline
bool
equal_xyz(const PointC3<K> &p, const PointC3<K> &q)
{
  return p.x() == q.x() && p.y() == q.y() && p.z() == q.z();
}

template < class K >
inline
Comparison_result
compare_xy(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().compare_xy_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_lexicographically_xy(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().compare_xy_3_object()(p, q);
}

template < class K >
inline
bool
lexicographically_xy_smaller_or_equal(const PointC3<K> &p, 
				      const PointC3<K> &q)
{ 
  return compare_lexicographically_xy(p, q) != LARGER;
}

template < class K >
inline
bool
lexicographically_xy_smaller(const PointC3<K> &p, 
			     const PointC3<K> &q)
{ 
  return K().less_xy_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_xyz(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().compare_xyz_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_lexicographically_xyz(const PointC3<K> &p,
                              const PointC3<K> &q)
{
  return K().compare_xyz_3_object()(p, q);
}

template < class K >
bool
lexicographically_xyz_smaller_or_equal(const PointC3<K> &p,
                                       const PointC3<K> &q)
{
  return compare_lexicographically_xyz(p, q) != LARGER;
}

template < class K >
inline
bool
lexicographically_xyz_smaller(const PointC3<K> &p,
                              const PointC3<K> &q)
{
  return K().less_xyz_3_object()(p, q);
}

template < class K >
inline
bool
strict_dominance(const PointC3<K> &p,
		 const PointC3<K> &q)
{
  return strict_dominanceC3(p.x(), p.y(), p.z(),
			    q.x(), q.y(), q.z());
}

template < class K >
inline
bool
dominance(const PointC3<K> &p,
	  const PointC3<K> &q)
{
  return dominanceC3(p.x(), p.y(), p.z(),
		     q.x(), q.y(), q.z());
}

template < class K >
inline
bool
collinear(const PointC3<K> &p,
          const PointC3<K> &q,
          const PointC3<K> &r)
{
    return K().collinear_3_object()(p, q, r);
}

template < class K >
inline
Orientation
orientation(const PointC3<K> &p,
            const PointC3<K> &q,
            const PointC3<K> &r,
            const PointC3<K> &s)
{
  return K().orientation_3_object()(p, q, r, s);
}

template < class K >
inline
Angle
angle(const PointC3<K> &p,
      const PointC3<K> &q,
      const PointC3<K> &r)
{
  return K().angle_3_object()(p, q, r);
}

template < class K >
inline
bool
coplanar(const PointC3<K> &p,
         const PointC3<K> &q,
         const PointC3<K> &r,
         const PointC3<K> &s)
{
  return K().coplanar_3_object()(p, q, r, s);
}

template < class K >
inline
Orientation
coplanar_orientation(const PointC3<K> &p,
                     const PointC3<K> &q,
                     const PointC3<K> &r,
                     const PointC3<K> &s)
{
  return K().coplanar_orientation_3_object()(p, q, r, s);
}

template < class K >
inline
Orientation
coplanar_orientation(const PointC3<K> &p,
                     const PointC3<K> &q,
                     const PointC3<K> &r)
{
  return K().coplanar_orientation_3_object()(p, q, r);
}

template < class K >
inline
Bounded_side
coplanar_side_of_bounded_circle(const PointC3<K> &p,
                                const PointC3<K> &q,
                                const PointC3<K> &r,
                                const PointC3<K> &t)
{
  return K().coplanar_side_of_bounded_circle_3_object()(p, q, r, t);
}

template < class K>
inline
bool
are_positive_oriented(const PointC3<K>& p,
                      const PointC3<K>& q,
                      const PointC3<K>& r,
                      const PointC3<K>& s)
{
  return orientation(p, q, r, s) == POSITIVE;
}

template < class K>
inline
bool
are_negative_oriented(const PointC3<K>& p,
                      const PointC3<K>& q,
                      const PointC3<K>& r,
                      const PointC3<K>& s)
{
  return orientation(p, q, r, s) == NEGATIVE;
}

template < class K >
inline
bool
are_ordered_along_line(const PointC3<K> &p,
                       const PointC3<K> &q,
                       const PointC3<K> &r)
{
  return K().are_ordered_along_line_3_object()(p, q, r);
}

template < class K >
inline
bool
collinear_are_ordered_along_line(const PointC3<K> &p,
                                 const PointC3<K> &q,
                                 const PointC3<K> &r)
{
  return K().collinear_are_ordered_along_line_3_object()(p, q, r);
}

template < class K >
inline
bool
are_strictly_ordered_along_line(const PointC3<K> &p,
                                const PointC3<K> &q,
                                const PointC3<K> &r)
{
  return K().are_strictly_ordered_along_line_3_object()(p, q, r);
}

template < class K >
inline
bool
collinear_are_strictly_ordered_along_line(const PointC3<K> &p,
                                          const PointC3<K> &q,
                                          const PointC3<K> &r)
{
  return K().collinear_are_strictly_ordered_along_line_3_object()(p, q, r);
}

template <class K >
inline
Oriented_side
side_of_oriented_sphere(const PointC3<K> &p,
                        const PointC3<K> &q,
                        const PointC3<K> &r,
                        const PointC3<K> &s,
                        const PointC3<K> &test)
{
  return K().side_of_oriented_sphere_3_object()(p, q, r, s, test);
}

template <class K >
inline
Bounded_side
side_of_bounded_sphere(const PointC3<K> &p,
                       const PointC3<K> &q,
                       const PointC3<K> &test)
{
  return K().side_of_bounded_sphere_3_object()(p, q, test);
}

template <class K >
inline
Bounded_side
side_of_bounded_sphere(const PointC3<K> &p,
                       const PointC3<K> &q,
                       const PointC3<K> &r,
                       const PointC3<K> &test)
{
  return K().side_of_bounded_sphere_3_object()(p, q, r, test);
}

template <class K >
inline
Bounded_side
side_of_bounded_sphere(const PointC3<K> &p,
                       const PointC3<K> &q,
                       const PointC3<K> &r,
                       const PointC3<K> &s,
                       const PointC3<K> &test)
{
  return K().side_of_bounded_sphere_3_object()(p, q, r, s, test);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_POINTS_3_H

// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_geom_traits_3.h
// revision      : $Revision$
// 
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================
//
// geometric traits for a <=3 D triangulation
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_GEOM_TRAITS_3_H
#define CGAL_TRIANGULATION_GEOM_TRAITS_3_H

#include <CGAL/basic.h>

#include <CGAL/Point_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class Repres >
class Triangulation_geom_traits_3 : public Repres
{
public:
  typedef Repres Rep;

  typedef typename Rep::Point_3        Point_3;
  typedef typename Rep::Segment_3      Segment_3;
  typedef typename Rep::Triangle_3     Triangle_3;
  typedef typename Rep::Tetrahedron_3  Tetrahedron_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3                      Point; 

  typedef typename Rep::Compare_x_3                Compare_x_3;
  typedef typename Rep::Compare_y_3                Compare_y_3;
  typedef typename Rep::Compare_z_3                Compare_z_3;
  typedef typename Rep::Equal_3                    Equal_3;
  typedef typename Rep::Collinear_3                Collinear_3;
  typedef typename Rep::Orientation_3              Orientation_3;
  typedef typename Rep::Coplanar_orientation_3     Coplanar_orientation_3;
  typedef typename Rep::Side_of_oriented_sphere_3  Side_of_oriented_sphere_3;
  typedef typename Rep::Coplanar_side_of_bounded_circle_3
                                          Coplanar_side_of_bounded_circle_3;

  typedef typename Rep::Construct_segment_3        Construct_segment_3;
  typedef typename Rep::Construct_triangle_3       Construct_triangle_3;
  typedef typename Rep::Construct_tetrahedron_3    Construct_tetrahedron_3;

  // For the hierarchy :
  typedef typename Rep::Less_distance_to_point_3 Less_distance_to_point_3;
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_GEOM_TRAITS_3_H

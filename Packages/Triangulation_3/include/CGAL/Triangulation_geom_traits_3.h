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

#include <CGAL/Cartesian.h>

#include <CGAL/Point_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/triangulation_assertions.h>

#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

// The classes prefixed with Local_ ae only needed because the
// kernel does not provide them yet.

template < class FT >
CGAL_KERNEL_LARGE_INLINE
Oriented_side
local_coplanar_side_of_oriented_circleC3(
	                           const FT &px, const FT &py, const FT &pz,
                                   const FT &qx, const FT &qy, const FT &qz,
                                   const FT &rx, const FT &ry, const FT &rz,
                                   const FT &tx, const FT &ty, const FT &tz,
                                   const FT &vx, const FT &vy, const FT &vz)
{
    // The approach is to compute side_of_oriented_sphere(p,q,r,t+v,t),
    // and remark that this expression simplifies internally.
  FT ptx = px - tx;
  FT pty = py - ty;
  FT ptz = pz - tz;
  FT pt2 = CGAL_NTS square(ptx) + CGAL_NTS square(pty) + CGAL_NTS square(ptz);
  FT qtx = qx - tx;
  FT qty = qy - ty;
  FT qtz = qz - tz;
  FT qt2 = CGAL_NTS square(qtx) + CGAL_NTS square(qty) + CGAL_NTS square(qtz);
  FT rtx = rx - tx;
  FT rty = ry - ty;
  FT rtz = rz - tz;
  FT rt2 = CGAL_NTS square(rtx) + CGAL_NTS square(rty) + CGAL_NTS square(rtz);
  FT v2 = CGAL_NTS square(vx) + CGAL_NTS square(vy) + CGAL_NTS square(vz);
  return Oriented_side(sign_of_determinant4x4(ptx,pty,ptz,pt2,
                                              rtx,rty,rtz,rt2,
                                              qtx,qty,qtz,qt2,
                                              vx,vy,vz,v2));
}

class Local_Coplanar_side_of_oriented_circle
{
  public:
    typedef Oriented_side  result_type;

    template <class T1, class T2>
    Oriented_side
    operator()(const T1& p, const T1& q, const T1& r, const T1& t,
	       const T2& v) const
    { 
      //return coplanar_side_of_oriented_circle(p,q,r,t, v); 
      return local_coplanar_side_of_oriented_circleC3(p.x(), p.y(), p.z(),
						      q.x(), q.y(), q.z(),
						      r.x(), r.y(), r.z(),
						      t.x(), t.y(), t.z(),
						      v.x(), v.y(), v.z());
    }
};

template < class Repres >
class Triangulation_geom_traits_3 
{
public:
  typedef Repres Rep;

  typedef typename Rep::Point_3        Point_3;
  typedef typename Rep::Segment_3      Segment_3;
  typedef typename Rep::Triangle_3     Triangle_3;
  typedef typename Rep::Tetrahedron_3  Tetrahedron_3;
  typedef typename Rep::Vector_3       Vector_3;

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
  typedef typename Rep::Construct_cross_product_vector_3
                                            Construct_cross_product_vector_3;

  // Uncomment the next line as soon as Kernels have this function
  // typedef typename Rep::Side_of_oriented_circle_3
  //                  Coplanar_side_of_oriented_circle_3;
  typedef Local_Coplanar_side_of_oriented_circle
                                          Coplanar_side_of_oriented_circle_3;

  typedef typename Rep::Construct_segment_3        Construct_segment_3;
  typedef typename Rep::Construct_triangle_3       Construct_triangle_3;
  typedef typename Rep::Construct_tetrahedron_3    Construct_tetrahedron_3;

  // And for the hierarchy :
  typedef typename Rep::Less_distance_to_point_3 Less_distance_to_point_3;

  Compare_x_3
  compare_x_3_object() const { 
    return Compare_x_3();
  }

  Compare_y_3
  compare_y_3_object() const { 
    return Compare_y_3();
  }

  Compare_z_3
  compare_z_3_object() const {
    return Compare_z_3();
  }
  
  Equal_3
  equal_3_object() const {
    return Equal_3();
  }

  Construct_cross_product_vector_3
  construct_cross_product_vector_3_object() const { 
    return Construct_cross_product_vector_3();
  }

  Collinear_3
  collinear_3_object() const {
    return Collinear_3();
  }

  Orientation_3
  orientation_3_object() const {
    return Orientation_3();
  }

  Coplanar_orientation_3
  coplanar_orientation_3_object() const {
    return Coplanar_orientation_3();
  }

  Coplanar_side_of_oriented_circle_3
  coplanar_side_of_oriented_circle_3_object() const {
    return Coplanar_side_of_oriented_circle_3();
  }

  Side_of_oriented_sphere_3
  side_of_oriented_sphere_3_object() const {
    return Side_of_oriented_sphere_3();
  }
 
  Construct_segment_3  construct_segment_3_object() const {
    return Construct_segment_3();
  }

  Construct_triangle_3  construct_triangle_3_object() const {
    return Construct_triangle_3();
  }

  Construct_tetrahedron_3  construct_tetrahedron_3_object() const {
    return Construct_tetrahedron_3();
  }

  Less_distance_to_point_3
  less_distance_to_point_3_object(const Point_3 &p) const
  {
      return Less_distance_to_point_3(p);
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_GEOM_TRAITS_3_H

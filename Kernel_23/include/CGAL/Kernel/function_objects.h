// Copyright (c) 1999,2002,2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany)
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stefan Schirra, Sylvain Pion

#ifndef CGAL_KERNEL_FUNCTION_OBJECTS_H
#define CGAL_KERNEL_FUNCTION_OBJECTS_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/intersection_2.h>
#include <CGAL/intersection_3.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL {

namespace CommonKernelFunctors {

  template <typename K>
  class Are_ordered_along_line_2
  {
    typedef typename K::Point_2     Point_2;
    typedef typename K::Collinear_2 Collinear_2;
    typedef typename K::Collinear_are_ordered_along_line_2
    Collinear_are_ordered_along_line_2;

    Collinear_2 c;
    Collinear_are_ordered_along_line_2 cao;
  public:
    typedef typename K::Boolean     result_type;

    Are_ordered_along_line_2() {}
    Are_ordered_along_line_2(const Collinear_2& c_,
			     const Collinear_are_ordered_along_line_2& cao_)
      : c(c_), cao(cao_)
    {}

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return c(p, q, r) && cao(p, q, r); }
  };

  template <typename K>
  class Are_ordered_along_line_3
  {
    typedef typename K::Point_3     Point_3;
    typedef typename K::Collinear_3 Collinear_3;
    typedef typename K::Collinear_are_ordered_along_line_3
    Collinear_are_ordered_along_line_3;

    Collinear_3 c;
    Collinear_are_ordered_along_line_3 cao;
  public:
    typedef typename K::Boolean     result_type;

    Are_ordered_along_line_3() {}
    Are_ordered_along_line_3(const Collinear_3& c_,
			     const Collinear_are_ordered_along_line_3& cao_)
      : c(c_), cao(cao_)
    {}

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return c(p, q, r) && cao(p, q, r); }
  };

  template <typename K>
  class Are_strictly_ordered_along_line_2
  {
    typedef typename K::Point_2     Point_2;
    typedef typename K::Collinear_2 Collinear_2;
    typedef typename K::Collinear_are_strictly_ordered_along_line_2
    Collinear_are_strictly_ordered_along_line_2;

    Collinear_2 c;
    Collinear_are_strictly_ordered_along_line_2 cao;
  public:
    typedef typename K::Boolean     result_type;

    Are_strictly_ordered_along_line_2() {}
    Are_strictly_ordered_along_line_2(
				      const Collinear_2& c_,
				      const Collinear_are_strictly_ordered_along_line_2& cao_)
      : c(c_), cao(cao_)
    {}

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return c(p, q, r) && cao(p, q, r); }
  };

  template <typename K>
  class Are_strictly_ordered_along_line_3
  {
    typedef typename K::Point_3     Point_3;
    typedef typename K::Collinear_3 Collinear_3;
    typedef typename K::Collinear_are_strictly_ordered_along_line_3
    Collinear_are_strictly_ordered_along_line_3;

    Collinear_3 c;
    Collinear_are_strictly_ordered_along_line_3 cao;
  public:
    typedef typename K::Boolean     result_type;

    Are_strictly_ordered_along_line_3() {}
    Are_strictly_ordered_along_line_3(
				      const Collinear_3& c_,
				      const Collinear_are_strictly_ordered_along_line_3& cao_)
      : c(c_), cao(cao_)
    {}

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return c(p, q, r) && cao(p, q, r); }
  };

  template <typename K>
  class Assign_2
  {
    typedef typename K::Object_2  Object_2;
  public:
    //typedef typename K::Boolean   result_type;
    typedef bool                  result_type;

    template <class T>
    result_type
    operator()(T& t, const Object_2& o) const
    { return assign(t, o); }
  };

  template <typename K>
  class Assign_3
  {
    typedef typename K::Object_3        Object_3;
  public:
    //typedef typename K::Boolean         result_type;
    typedef bool                        result_type;

    template <class T>
    result_type
    operator()(T& t, const Object_3& o) const
    { return assign(t, o); }
  };

  template <typename K>
  class Compare_dihedral_angle_3
  {
    typedef typename K::Point_3            Point_3;
    typedef typename K::Vector_3           Vector_3;
    typedef typename K::FT                 FT;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()(const Point_3& a1, const Point_3& b1, 
               const Point_3& c1, const Point_3& d1, 
               const Point_3& a2, const Point_3& b2, 
               const Point_3& c2, const Point_3& d2) const
    {
      const Vector_3 ab1 = b1 - a1;
      const Vector_3 ac1 = c1 - a1;
      const Vector_3 ad1 = d1 - a1;

      const Vector_3 ab2 = b2 - a2;
      const Vector_3 ac2 = c2 - a2;
      const Vector_3 ad2 = d2 - a2;
      return this->operator()(ab1, ac1, ad1, ab2, ac2, ad2);
    }

    result_type
    operator()(const Point_3& a1, const Point_3& b1, 
               const Point_3& c1, const Point_3& d1, 
               const FT& cosine) const
    {
      const Vector_3 ab1 = b1 - a1;
      const Vector_3 ac1 = c1 - a1;
      const Vector_3 ad1 = d1 - a1;

      return this->operator()(ab1, ac1, ad1, cosine);
    }

    result_type
    operator()(const Vector_3& ab1, const Vector_3& ac1, const Vector_3& ad1,
               const FT& cosine)
      const
    {
      typedef typename K::FT                                 FT;
      typedef typename K::Construct_cross_product_vector_3   Cross_product;
      Cross_product xproduct = K().construct_cross_product_vector_3_object();

      const Vector_3 abac1 = xproduct(ab1, ac1);
      const Vector_3 abad1 = xproduct(ab1, ad1);
      const FT sc_prod_1 = abac1 * abad1;

      CGAL_kernel_assertion_msg( abac1 != NULL_VECTOR,
                                 "ab1 and ac1 are collinear" );
      CGAL_kernel_assertion_msg( abad1 != NULL_VECTOR,
                                 "ab1 and ad1 are collinear" );

      if(sc_prod_1 >= 0 ) {
        if(cosine >= 0) {
          // the two cosine are >= 0, cosine is decreasing on [0,1]
          return compare(CGAL::square(cosine)*
                         abac1.squared_length()*abad1.squared_length(),
                         CGAL::square(sc_prod_1));
        }
        else {
          return SMALLER;
        }
      }
      else {
        if(cosine < 0) {
          // the two cosine are < 0, cosine is increasing on [-1,0]
          return compare(CGAL::square(sc_prod_1),
                         CGAL::square(cosine)*
                         abac1.squared_length()*abad1.squared_length());
        }
        else
          return LARGER;
        }
    }

    result_type
    operator()(const Vector_3& ab1, const Vector_3& ac1, const Vector_3& ad1,
               const Vector_3& ab2, const Vector_3& ac2, const Vector_3& ad2)
      const
    {
      typedef typename K::FT                                 FT;
      typedef typename K::Construct_cross_product_vector_3   Cross_product;
      Cross_product xproduct = K().construct_cross_product_vector_3_object();

      const Vector_3 abac1 = xproduct(ab1, ac1);
      const Vector_3 abad1 = xproduct(ab1, ad1);
      const FT sc_prod_1 = abac1 * abad1;

      const Vector_3 abac2 = xproduct(ab2, ac2);
      const Vector_3 abad2 = xproduct(ab2, ad2);
      const FT sc_prod_2 = abac2 * abad2;

      CGAL_kernel_assertion_msg( abac1 != NULL_VECTOR,
                                 "ab1 and ac1 are collinear" );
      CGAL_kernel_assertion_msg( abad1 != NULL_VECTOR,
                                 "ab1 and ad1 are collinear" );
      CGAL_kernel_assertion_msg( abac2 != NULL_VECTOR,
                                 "ab2 and ac2 are collinear" );
      CGAL_kernel_assertion_msg( abad2 != NULL_VECTOR,
                                 "ab2 and ad2 are collinear" );

      if(sc_prod_1 >= 0 ) {
        if(sc_prod_2 >= 0) {
          // the two cosine are >= 0, cosine is decreasing on [0,1]
          return compare(CGAL::square(sc_prod_2)*
                         abac1.squared_length()*abad1.squared_length(),
                         CGAL::square(sc_prod_1)*
                         abac2.squared_length()*abad2.squared_length());
        }
        else {
          return SMALLER;
        }
      }
      else {
        if(sc_prod_2 < 0) {
          // the two cosine are < 0, cosine is increasing on [-1,0]
          return compare(CGAL::square(sc_prod_1)*
                         abac2.squared_length()*abad2.squared_length(),
                         CGAL::square(sc_prod_2)*
                         abac1.squared_length()*abad1.squared_length());
        }
        else
          return LARGER;
        }
    }
  };

  template <typename K>
  class Compare_squared_distance_2
  {
    typedef typename K::FT                 FT;
  public:
    typedef typename K::Comparison_result  result_type;

    template <class T1, class T2>
    result_type
    operator()(const T1& p, const T2& q, const FT& d2) const
    {
      return CGAL_NTS compare(squared_distance(p, q), d2);
    }

    template <class T1, class T2, class T3, class T4>
    result_type
    operator()(const T1& p, const T2& q, const T3& r, const T4& s) const
    {
      return CGAL_NTS compare(squared_distance(p, q), squared_distance(r, s));
    }
  };

  template <typename K>
  class Compare_squared_distance_3
  {
    typedef typename K::FT                 FT;
  public:
    typedef typename K::Comparison_result  result_type;

    template <class T1, class T2>
    result_type
    operator()(const T1& p, const T2& q, const FT& d2) const
    {
      return CGAL_NTS compare(squared_distance(p, q), d2);
    }

    template <class T1, class T2, class T3, class T4>
    result_type
    operator()(const T1& p, const T2& q, const T3& r, const T4& s) const
    {
      return CGAL_NTS compare(squared_distance(p, q), squared_distance(r, s));
    }
  };

  template <typename K>
  class Compute_area_3
  {
    typedef typename K::FT                FT;
    typedef typename K::Point_3           Point_3;
    typedef typename K::Triangle_3        Triangle_3;
  public:
    typedef FT               result_type;

    FT
    operator()( const Triangle_3& t ) const
    {
	return CGAL_NTS sqrt(K().compute_squared_area_3_object()(t));
    }

    FT
    operator()( const Point_3& p, const Point_3& q, const Point_3& r ) const
    {
	return CGAL_NTS sqrt(K().compute_squared_area_3_object()(p, q, r));
    }
  };

  template <typename K>
  class Compute_squared_distance_2
  {
    typedef typename K::FT   FT;
  public:
    typedef FT               result_type;

    // There are 25 combinaisons, we use a template.
    template <class T1, class T2>
    FT
    operator()( const T1& t1, const T2& t2) const
    { return internal::squared_distance(t1, t2, K()); }
  };

  template <typename K>
  class Compute_squared_distance_3
  {
    typedef typename K::FT        FT;
    typedef typename K::Point_3   Point_3;
  public:
    typedef FT               result_type;

    // There are 25 combinaisons, we use a template.
    template <class T1, class T2>
    FT
    operator()( const T1& t1, const T2& t2) const
    { return internal::squared_distance(t1, t2, K()); }

    FT
    operator()( const Point_3& pt1, const Point_3& pt2) const
    {
      typedef typename K::Vector_3 Vector_3;
      Vector_3 vec = pt2 - pt1;
      return vec*vec;
    }
  };

  template <typename K>
  class Compute_squared_length_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Segment_2   Segment_2;
    typedef typename K::Vector_2    Vector_2;
  public:
    typedef FT               result_type;

    FT
    operator()( const Vector_2& v) const
    { return CGAL_NTS square(K().compute_x_2_object()(v)) +
             CGAL_NTS square(K().compute_y_2_object()(v));}

    FT
    operator()( const Segment_2& s) const
    { return K().compute_squared_distance_2_object()(s.source(), s.target()); }
  };

  template <typename K>
  class Compute_squared_length_3
  {
    typedef typename K::FT          FT;
    typedef typename K::Segment_3   Segment_3;
    typedef typename K::Vector_3    Vector_3;
  public:
    typedef FT               result_type;

    FT
    operator()( const Vector_3& v) const
    { return v.rep().squared_length(); }

    FT
    operator()( const Segment_3& s) const
    { return s.squared_length(); }
  };

  template <typename K>
  class Compute_a_2
  {
    typedef typename K::RT             RT;
    typedef typename K::Line_2         Line_2;

  public:
    typedef RT               result_type;

    RT
    operator()(const Line_2& l) const
    {
      return l.rep().a();
    }
  };

  template <typename K>
  class Compute_a_3
  {
    typedef typename K::RT             RT;
    typedef typename K::Plane_3        Plane_3;

  public:
    typedef RT               result_type;

    RT
    operator()(const Plane_3& l) const
    {
      return l.rep().a();
    }
  };


  template <typename K>
  class Compute_b_2
  {
    typedef typename K::RT             RT;
    typedef typename K::Line_2         Line_2;

  public:
    typedef RT               result_type;

    RT
    operator()(const Line_2& l) const
    {
      return l.rep().b();
    }
  };

  template <typename K>
  class Compute_b_3
  {
    typedef typename K::RT             RT;
    typedef typename K::Plane_3        Plane_3;

  public:
    typedef RT               result_type;

    RT
    operator()(const Plane_3& l) const
    {
      return l.rep().b();
    }
  };


  template <typename K>
  class Compute_c_2
  {
    typedef typename K::RT             RT;
    typedef typename K::Line_2         Line_2;

  public:
    typedef RT               result_type;

    RT
    operator()(const Line_2& l) const
    {
      return l.rep().c();
    }
  };

  template <typename K>
  class Compute_c_3
  {
    typedef typename K::RT             RT;
    typedef typename K::Plane_3        Plane_3;

  public:
    typedef RT               result_type;

    RT
    operator()(const Plane_3& l) const
    {
      return l.rep().c();
    }
  };

  template <typename K>
  class Compute_d_3
  {
    typedef typename K::RT             RT;
    typedef typename K::Plane_3        Plane_3;

  public:
    typedef RT               result_type;

    RT
    operator()(const Plane_3& l) const
    {
      return l.rep().d();
    }
  };

  template <typename K>
  class Compute_x_at_y_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Line_2         Line_2;

  public:
    typedef FT               result_type;

    FT
    operator()(const Line_2& l, const FT& y) const
    {
      CGAL_kernel_precondition( ! l.is_degenerate() );
      return (FT(-l.b())*y - FT(l.c()) )/FT(l.a());
    }
  };

  template <typename K>
  class Compute_y_at_x_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Line_2         Line_2;

  public:
    typedef FT               result_type;

    FT
    operator()(const Line_2& l, const FT& x) const
    {
      CGAL_kernel_precondition_msg( ! l.is_vertical(),
		    "Compute_y_at_x(FT x) is undefined for vertical line");
      return (FT(-l.a())*x - FT(l.c()) )/FT(l.b());
    }
  };

  template <typename K>
  class Compute_xmin_2
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_rectangle_2& r) const
    {
      return (r.min)().x();
    }
  };

  template <typename K>
  class Compute_xmin_3
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_cuboid_3    Iso_cuboid_3;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_cuboid_3& r) const
    {
      return (r.min)().x();
    }
  };

  template <typename K>
  class Compute_xmax_2
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_rectangle_2& r) const
    {
      return (r.max)().x();
    }
  };

  template <typename K>
  class Compute_xmax_3
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_cuboid_3    Iso_cuboid_3;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_cuboid_3& r) const
    {
      return (r.max)().x();
    }
  };

  template <typename K>
  class Compute_ymin_2
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_rectangle_2& r) const
    {
      return (r.min)().y();
    }
  };

  template <typename K>
  class Compute_ymin_3
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_cuboid_3    Iso_cuboid_3;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_cuboid_3& r) const
    {
      return (r.min)().y();
    }
  };

  template <typename K>
  class Compute_ymax_2
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_rectangle_2& r) const
    {
      return (r.max)().y();
    }
  };

  template <typename K>
  class Compute_ymax_3
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_cuboid_3    Iso_cuboid_3;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_cuboid_3& r) const
    {
      return (r.max)().y();
    }
  };

  template <typename K>
  class Compute_zmin_3
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_cuboid_3    Iso_cuboid_3;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_cuboid_3& r) const
    {
      return (r.min)().z();
    }
  };

  template <typename K>
  class Compute_zmax_3
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_cuboid_3    Iso_cuboid_3;
    typedef FT                          Cartesian_coordinate_type;
    //typedef typename K::Cartesian_coordinate_type  Cartesian_coordinate_type;

  public:
    typedef FT               result_type;

    Cartesian_coordinate_type
    operator()(const Iso_cuboid_3& r) const
    {
      return (r.max)().z();
    }
  };

  template <typename K>
  class Construct_center_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Circle_2  Circle_2;
  public:
    typedef const Point_2&         result_type;

    result_type
    operator()(const Circle_2& c) const
    { return c.rep().center(); }
  };

  template <typename K>
  class Construct_center_3
  {
    typedef typename K::Point_3   Point_3;
    typedef typename K::Sphere_3  Sphere_3;
    typedef typename K::Circle_3  Circle_3;
  public:
    typedef const Point_3&          result_type;

    result_type
    operator()(const Sphere_3& s) const
    { return s.rep().center(); }

    result_type
    operator()(const Circle_3& c) const
    { return c.rep().center(); }

  };

  template <typename K>
  class Construct_circle_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Circle_2    Circle_2;
    typedef typename Circle_2::Rep  Rep;
  public:
    typedef Circle_2         result_type;

    Rep // Circle_2
    operator()( Return_base_tag,
                const Point_2& center, const FT& squared_radius,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return Rep(center, squared_radius, orientation); }

    Rep // Circle_2
    operator()( Return_base_tag,
                const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      typename K::Orientation_2 orientation;
      typename K::Compute_squared_distance_2 squared_distance;
      typename K::Construct_circumcenter_2 circumcenter;
      typename K::Orientation orient = orientation(p, q, r);
      CGAL_kernel_precondition( orient != COLLINEAR);

      Point_2 center = circumcenter(p, q, r);

      return Rep(center, squared_distance(p, center), orient);
    }

    Rep // Circle_2
    operator()( Return_base_tag,
                const Point_2& p, const Point_2& q,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {
      CGAL_kernel_precondition( orientation != COLLINEAR);

      typename K::Compute_squared_distance_2 squared_distance;
      typename K::Construct_midpoint_2 midpoint;
      if (p != q) {
        Point_2 center = midpoint(p, q);
        return Rep(center, squared_distance(p, center), orientation);
      } else
        return Rep(p, FT(0), orientation);
    }

    Rep // Circle_2
    operator()( Return_base_tag,
                const Point_2& p, const Point_2& q,
	        const FT& bulge) const
    {
     
      typename K::Compute_squared_distance_2 squared_distance;
      const FT sqr_bulge = CGAL::square(bulge);
      const FT common = (FT(1) - sqr_bulge) / (FT(4)*bulge);
      const FT x_coord = (p.x() + q.x())/FT(2)
	                 + common*(p.y() - q.y());
      const FT y_coord = (p.y() + q.y())/FT(2)
                          + common*(q.x() - p.x());
      
      const FT sqr_rad = squared_distance(p, q) 
	                 * (FT(1)/sqr_bulge + FT(2) + sqr_bulge) / FT(16); 

      return Rep(Point_2(x_coord, y_coord), sqr_rad); 
    }


    Rep // Circle_2
    operator()( Return_base_tag, const Point_2& center,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {
      CGAL_kernel_precondition( orientation != COLLINEAR );

      return Rep(center, FT(0), orientation);
    }


    Circle_2
    operator()( const Point_2& center, const FT& squared_radius,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {
      return this->operator()(Return_base_tag(),
                              center, squared_radius, orientation);
    }

    Circle_2
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      return this->operator()(Return_base_tag(), p, q, r);
    }

    Circle_2
    operator()( const Point_2& p, const Point_2& q,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {
      return this->operator()(Return_base_tag(), p, q, orientation);
    }

    Circle_2
    operator()( const Point_2& p, const Point_2& q,
	        const FT& bulge) const
    {
      return this->operator()(Return_base_tag(), p, q, bulge);
    }

    Circle_2
    operator()( const Point_2& center,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {
      return this->operator()(Return_base_tag(), center, orientation);
    }
  };

  template < typename K >
  class Construct_circle_3
  {
    typedef typename K::FT           FT;
    typedef typename K::Point_3      Point_3;
    typedef typename K::Plane_3      Plane_3;
    typedef typename K::Sphere_3     Sphere_3;
    typedef typename K::Circle_3     Circle_3;
    typedef typename K::Vector_3     Vector_3;
    typedef typename K::Direction_3  Direction_3;
    typedef typename Circle_3::Rep    Rep;

  public:
    typedef Circle_3                  result_type;

    Rep
    operator() (Return_base_tag, const Point_3& p,
                const FT& sr, const Plane_3& plane) const
    { return Rep(p, sr, plane); }

    Rep
    operator() (Return_base_tag, const Point_3& p,
                const FT& sr, const Vector_3& v) const
    { return Rep(p, sr, v); }

    Rep
    operator() (Return_base_tag, const Point_3& p,
                const FT& sr, const Direction_3& d) const
    { return Rep(p, sr, d); }

    Rep
    operator() (Return_base_tag, const Sphere_3& s1,
                const Sphere_3& s2) const
    { return Rep(s1, s2); }

    Rep
    operator() (Return_base_tag, const Plane_3& p,
                const Sphere_3& s) const
    { return Rep(p, s); }

    Rep
    operator() (Return_base_tag, const Plane_3& p,
                const Sphere_3& s, int a) const
    { return Rep(p, s, a); }

    Rep
    operator() (Return_base_tag, const Point_3& p1,
                const Point_3& p2, const Point_3& p3) const
    { return Rep(p1, p2, p3); }

    Circle_3
    operator()(const Point_3& p, const FT& sr,
               const Plane_3& plane) const
    { return this->operator()(Return_base_tag(), p, sr, plane); }

    Circle_3
    operator() (const Point_3& p, const FT& sr,
                const Vector_3& v) const
    { return this->operator()(Return_base_tag(), p, sr, v); }

    Circle_3
    operator() (const Point_3& p, const FT& sr,
                const Direction_3& d) const
    { return this->operator()(Return_base_tag(), p, sr, d); }

    Circle_3
    operator() (const Sphere_3& s1, const Sphere_3& s2) const
    { return this->operator()(Return_base_tag(), s1, s2); }

    Circle_3
    operator() (const Plane_3& p, const Sphere_3& s) const
    { return this->operator()(Return_base_tag(), p, s); }

    Circle_3
    operator() (const Sphere_3& s, const Plane_3& p) const
    { return this->operator()(Return_base_tag(), p, s); }

    Circle_3
    operator() (const Plane_3& p, const Sphere_3& s, int a) const
    { return this->operator()(Return_base_tag(), p, s, a); }

    Circle_3
    operator() (const Sphere_3& s, const Plane_3& p, int a) const
    { return this->operator()(Return_base_tag(), p, s, a); }

    Circle_3
    operator()(	const Point_3& p1, const Point_3& p2, const Point_3& p3) const
    { return this->operator()(Return_base_tag(), p1, p2, p3); }
  };

  template <typename K>
  class Construct_iso_cuboid_3
  {
    typedef typename K::RT            RT;
    typedef typename K::Point_3       Point_3;
    typedef typename K::Iso_cuboid_3  Iso_cuboid_3;
    typedef typename Iso_cuboid_3::Rep  Rep;
  public:
    typedef Iso_cuboid_3      result_type;

    Rep // Iso_cuboid_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q, int) const
    { return Rep(p, q, 0); }

    Rep // Iso_cuboid_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q) const
    { return Rep(p, q); }

    Rep // Iso_cuboid_3
    operator()(Return_base_tag, const Point_3 &left,   const Point_3 &right,
               const Point_3 &bottom, const Point_3 &top,
               const Point_3 &far_,   const Point_3 &close) const
    { return Rep(left, right, bottom, top, far_, close); }

    Rep // Iso_cuboid_3
    operator()(Return_base_tag, const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz,
               const RT& hw) const
    { return Rep(min_hx, min_hy, min_hz, max_hx, max_hy, max_hz, hw); }

    Rep // Iso_cuboid_3
    operator()(Return_base_tag, const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz) const
    { return Rep(min_hx, min_hy, min_hz, max_hx, max_hy, max_hz); }


    Iso_cuboid_3
    operator()(const Point_3& p, const Point_3& q, int) const
    { return this->operator()(Return_base_tag(), p, q, 0); }

    Iso_cuboid_3
    operator()(const Point_3& p, const Point_3& q) const
    { return this->operator()(Return_base_tag(), p, q); }

    Iso_cuboid_3
    operator()(const Point_3 &left,   const Point_3 &right,
               const Point_3 &bottom, const Point_3 &top,
               const Point_3 &far_,   const Point_3 &close) const
    { return this->operator()(Return_base_tag(), left, right, bottom, top, far_, close); }

    Iso_cuboid_3
    operator()(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz,
               const RT& hw) const
    { return this->operator()(Return_base_tag(), min_hx, min_hy, min_hz, max_hx, max_hy, max_hz, hw); }

    Iso_cuboid_3
    operator()(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz) const
    { return this->operator()(Return_base_tag(), min_hx, min_hy, min_hz, max_hx, max_hy, max_hz); }
  };

  template <typename K>
  class Construct_max_vertex_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Segment_2        Segment_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  public:
    typedef const Point_2&               result_type;

    result_type
    operator()(const Iso_rectangle_2& r) const
    { return (r.rep().max)(); }

    result_type
    operator()(const Segment_2& s) const
    { return (s.max)(); }
  };


  template <typename K>
  class Construct_min_vertex_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Segment_2        Segment_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  public:
    typedef const Point_2&               result_type;

    result_type
    operator()(const Iso_rectangle_2& r) const
    { return (r.rep().min)(); }

    result_type
    operator()(const Segment_2& s) const
    { return (s.min)(); }
  };



  template <typename K>
  class Construct_max_vertex_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
  public:
    typedef const Point_3&           result_type;

    result_type
    operator()(const Iso_cuboid_3& r) const
    { return (r.rep().max)(); }

    result_type
    operator()(const Segment_3& s) const
    { return (s.rep().max)(); }
  };

  template <typename K>
  class Construct_min_vertex_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
  public:
    typedef const Point_3&               result_type;

    result_type
    operator()(const Iso_cuboid_3& r) const
    { return (r.rep().min)(); }

    result_type
    operator()(const Segment_3& s) const
    { return (s.rep().min)(); }
  };

  template <typename K>
  class Construct_normal_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Vector_3         Vector_3;
  public:
    typedef Vector_3           result_type;

    Vector_3
    operator()(const Point_3& p,const Point_3& q, const Point_3& r) const
    { 
      CGAL_kernel_precondition(! K().collinear_3_object()(p,q,r) ); 
      Vector_3 res = CGAL::cross_product(q-p, r-p);
      return res; }
  };

  template <typename K>
  class Construct_object_2
  {
    typedef typename K::Object_2   Object_2;
  public:
    typedef Object_2         result_type;

    template <class Cls>
    Object_2
    operator()( const Cls& c) const
    { return make_object(c); }
  };

  template <typename K>
  class Construct_object_3
  {
    typedef typename K::Object_3   Object_3;
  public:
    typedef Object_3         result_type;

    template <class Cls>
    Object_3
    operator()( const Cls& c) const
    { return make_object(c); }
  };

  template <typename K>
  class Construct_opposite_circle_2
  {
    typedef typename K::Circle_2   Circle_2;
  public:
    typedef Circle_2         result_type;

    Circle_2
    operator()( const Circle_2& c) const
    { return c.opposite(); }
  };

  template <typename K>
  class Construct_opposite_direction_2
  {
    typedef typename K::Direction_2    Direction_2;
    typedef typename Direction_2::Rep  Rep;
  public:
    typedef Direction_2      result_type;

    Direction_2
    operator()( const Direction_2& d) const
    {  return Rep(-d.dx(), -d.dy()); }
  };

  template <typename K>
  class Construct_opposite_direction_3
  {
    typedef typename K::Direction_3    Direction_3;
    typedef typename Direction_3::Rep  Rep;
  public:
    typedef Direction_3      result_type;

    Direction_3
    operator()( const Direction_3& d) const
    {  return Rep(-d.dx(), -d.dy(), -d.dz()); }
  };

  template <typename K>
  class Construct_opposite_line_2
  {
    typedef typename K::Line_2   Line_2;
  public:
    typedef Line_2           result_type;

    Line_2
    operator()( const Line_2& l) const
    { return Line_2( -l.a(), -l.b(), -l.c()); }
  };

  template <typename K>
  class Construct_opposite_line_3
  {
    typedef typename K::Line_3   Line_3;
  public:
    typedef Line_3           result_type;

    Line_3
    operator()( const Line_3& l) const
    { return l.rep().opposite(); }
  };

  template <typename K>
  class Construct_opposite_plane_3
  {
    typedef typename K::Plane_3   Plane_3;
  public:
    typedef Plane_3          result_type;

    Plane_3
    operator()( const Plane_3& p) const
    { return p.rep().opposite(); }
  };

  template <typename K>
  class Construct_opposite_ray_2
  {
    typedef typename K::Ray_2   Ray_2;
  public:
    typedef Ray_2            result_type;

    Ray_2
    operator()( const Ray_2& r) const
    { return r.opposite(); }
  };

  template <typename K>
  class Construct_opposite_ray_3
  {
    typedef typename K::Ray_3   Ray_3;
  public:
    typedef Ray_3            result_type;

    Ray_3
    operator()( const Ray_3& r) const
    { return r.opposite(); }
  };

  template <typename K>
  class Construct_opposite_segment_2
  {
    typedef typename K::Segment_2  Segment_2;
  public:
    typedef Segment_2        result_type;

    Segment_2
    operator()( const Segment_2& s) const
    { return Segment_2(s.target(), s.source()); }
  };

  template <typename K>
  class Construct_opposite_segment_3
  {
    typedef typename K::Segment_3  Segment_3;
  public:
    typedef Segment_3        result_type;

    Segment_3
    operator()( const Segment_3& s) const
    { return s.rep().opposite(); }
  };

  template <typename K>
  class Construct_opposite_sphere_3
  {
    typedef typename K::Sphere_3   Sphere_3;
  public:
    typedef Sphere_3         result_type;

    Sphere_3
    operator()( const Sphere_3& s) const
    { return s.rep().opposite(); }
  };

  template <typename K>
  class Construct_opposite_triangle_2
  {
    typedef typename K::Triangle_2  Triangle_2;
  public:
    typedef Triangle_2       result_type;

    Triangle_2
    operator()( const Triangle_2& t) const
    { return Triangle_2(t.vertex(0), t.vertex(2), t.vertex(1));}
  };

  template <typename K>
  class Construct_perpendicular_line_3
  {
    typedef typename K::Line_3    Line_3;
    typedef typename K::Point_3   Point_3;
    typedef typename K::Plane_3   Plane_3;
  public:
    typedef Line_3           result_type;

    Line_3
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.rep().perpendicular_line(p); }
  };

  template <typename K>
  class Construct_perpendicular_plane_3
  {
    typedef typename K::Line_3    Line_3;
    typedef typename K::Point_3   Point_3;
    typedef typename K::Plane_3   Plane_3;
  public:
    typedef Plane_3          result_type;

    Plane_3
    operator()( const Line_3& l, const Point_3& p) const
    { return l.rep().perpendicular_plane(p); }
  };

  template <typename K>
  class Construct_plane_3
  {
    typedef typename K::RT           RT;
    typedef typename K::Point_3      Point_3;
    typedef typename K::Vector_3     Vector_3;
    typedef typename K::Direction_3  Direction_3;
    typedef typename K::Line_3       Line_3;
    typedef typename K::Ray_3        Ray_3;
    typedef typename K::Segment_3    Segment_3;
    typedef typename K::Plane_3      Plane_3;
    typedef typename K::Circle_3     Circle_3;
    typedef typename Plane_3::Rep    Rep;
  public:
    typedef Plane_3          result_type;

    Rep // Plane_3
    operator()(Return_base_tag, const RT& a, const RT& b, const RT& c, const RT& d) const
    { return Rep(a, b, c, d); }

    Rep // Plane_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q, const Point_3& r) const
    { return Rep(p, q, r); }

    Rep // Plane_3
    operator()(Return_base_tag, const Point_3& p, const Direction_3& d) const
    { return Rep(p, d); }

    Rep // Plane_3
    operator()(Return_base_tag, const Point_3& p, const Vector_3& v) const
    { return Rep(p, v); }

    Rep // Plane_3
    operator()(Return_base_tag, const Line_3& l, const Point_3& p) const
    { return Rep(l, p); }

    Rep // Plane_3
    operator()(Return_base_tag, const Ray_3& r, const Point_3& p) const
    { return Rep(r, p); }

    Rep // Plane_3
    operator()(Return_base_tag, const Segment_3& s, const Point_3& p) const
    { return Rep(s, p); }

    Rep // Plane_3
    operator()(Return_base_tag, const Circle_3 & c) const
    { return c.rep().supporting_plane(); }

    Plane_3
    operator()(const RT& a, const RT& b, const RT& c, const RT& d) const
    { return this->operator()(Return_base_tag(), a, b, c, d); }

    Plane_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return this->operator()(Return_base_tag(), p, q, r); }

    Plane_3
    operator()(const Point_3& p, const Direction_3& d) const
    { return this->operator()(Return_base_tag(), p, d); }

    Plane_3
    operator()(const Point_3& p, const Vector_3& v) const
    { return this->operator()(Return_base_tag(), p, v); }

    Plane_3
    operator()(const Line_3& l, const Point_3& p) const
    { return this->operator()(Return_base_tag(), l, p); }

    Plane_3
    operator()(const Ray_3& r, const Point_3& p) const
    { return this->operator()(Return_base_tag(), r, p); }

    Plane_3
    operator()(const Segment_3& s, const Point_3& p) const
    { return this->operator()(Return_base_tag(), s, p); }

    Plane_3
    operator()(const Circle_3 & c) const
    { return this->operator()(Return_base_tag(), c); }

  };

  template <typename K>
  class Construct_point_on_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Segment_2  Segment_2;
    typedef typename K::Line_2     Line_2;
    typedef typename K::Ray_2      Ray_2;
  public:
    typedef Point_2          result_type;

    Point_2
    operator()( const Line_2& l, int i) const
    { return l.point(i); }

    Point_2
    operator()( const Segment_2& s, int i) const
    { return s.point(i); }

    Point_2
    operator()( const Ray_2& r, int i) const
    { return r.point(i); }
  };

  template <typename K>
  class Construct_point_on_3
  {
    typedef typename K::Point_3    Point_3;
    typedef typename K::Segment_3  Segment_3;
    typedef typename K::Line_3     Line_3;
    typedef typename K::Ray_3      Ray_3;
    typedef typename K::Plane_3    Plane_3;
  public:
    typedef Point_3          result_type;

    Point_3
    operator()( const Line_3& l, int i) const
    { return l.rep().point(i); }

    Point_3
    operator()( const Segment_3& s, int i) const
    { return s.point(i); }

    Point_3
    operator()( const Ray_3& r, int i) const
    { return r.rep().point(i); }

    Point_3
    operator()( const Plane_3& p) const
    { return p.rep().point(); }
  };

  template <typename K>
  class Construct_projected_xy_point_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Point_3    Point_3;
    typedef typename K::Plane_3    Plane_3;
  public:
    typedef Point_2          result_type;

    Point_2
    operator()( const Plane_3& h, const Point_3& p) const
    {  return h.rep().to_2d(p); }
  };

  template <typename K>
  class Construct_ray_2
  {
    typedef typename K::Point_2      Point_2;
    typedef typename K::Vector_2     Vector_2;
    typedef typename K::Direction_2  Direction_2;
    typedef typename K::Line_2       Line_2;
    typedef typename K::Ray_2        Ray_2;
    typedef typename Ray_2::Rep   Rep;
  public:
    typedef Ray_2            result_type;

    Rep // Ray_2
    operator()(Return_base_tag, const Point_2& p, const Point_2& q) const
    {  return Rep(p, q); }

    Rep // Ray_2
    operator()(Return_base_tag, const Point_2& p, const Vector_2& v) const
    {  return Rep(p, K().construct_translated_point_2_object()(p,  v)); }

    Rep // Ray_2
    operator()(Return_base_tag, const Point_2& p, const Direction_2& d) const
    {  return Rep(p, K().construct_translated_point_2_object()(p, d.to_vector())); }

    Rep // Ray_2
    operator()(Return_base_tag, const Point_2& p, const Line_2& l) const
    {  return Rep(p, K().construct_translated_point_2_object()(p, l.to_vector())); }


    Ray_2
    operator()(const Point_2& p, const Point_2& q) const
    { return this->operator()(Return_base_tag(), p, q); }

    Ray_2
    operator()(const Point_2& p, const Vector_2& v) const
    { return this->operator()(Return_base_tag(), p, v); }

    Ray_2
    operator()(const Point_2& p, const Direction_2& d) const
    { return this->operator()(Return_base_tag(), p, d); }

    Ray_2
    operator()(const Point_2& p, const Line_2& l) const
    { return this->operator()(Return_base_tag(), p, l); }
  };

  template <typename K>
  class Construct_ray_3
  {
    typedef typename K::Point_3      Point_3;
    typedef typename K::Vector_3     Vector_3;
    typedef typename K::Direction_3  Direction_3;
    typedef typename K::Line_3       Line_3;
    typedef typename K::Ray_3        Ray_3;
    typedef typename Ray_3::Rep      Rep;
  public:
    typedef Ray_3            result_type;

    Rep // Ray_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q) const
    {  return Rep(p, q); }

    Rep // Ray_3
    operator()(Return_base_tag, const Point_3& p, const Vector_3& v) const
    {  return Rep(p, v); }

    Rep // Ray_3
    operator()(Return_base_tag, const Point_3& p, const Direction_3& d) const
    {  return Rep(p, d); }

    Rep // Ray_3
    operator()(Return_base_tag, const Point_3& p, const Line_3& l) const
    {  return Rep(p, l); }


    Ray_3
    operator()(const Point_3& p, const Point_3& q) const
    { return this->operator()(Return_base_tag(), p, q); }

    Ray_3
    operator()(const Point_3& p, const Vector_3& v) const
    { return this->operator()(Return_base_tag(), p, v); }

    Ray_3
    operator()(const Point_3& p, const Direction_3& d) const
    { return this->operator()(Return_base_tag(), p, d); }

    Ray_3
    operator()(const Point_3& p, const Line_3& l) const
    { return this->operator()(Return_base_tag(), p, l); }
  };

  template <typename K>
  class Construct_segment_2
  {
    typedef typename K::Segment_2  Segment_2;
    typedef typename Segment_2::Rep  Rep;
    typedef typename K::Point_2    Point_2;
  public:
    typedef Segment_2        result_type;

    Rep // Segment_2
    operator()(Return_base_tag, const Point_2& p, const Point_2& q) const
    {  return Rep(p, q); }

    Segment_2
    operator()( const Point_2& p, const Point_2& q) const
    { return this->operator()(Return_base_tag(), p, q); }
  };

  template <typename K>
  class Construct_segment_3
  {
    typedef typename K::Segment_3  Segment_3;
    typedef typename K::Point_3    Point_3;
    typedef typename Segment_3::Rep  Rep;
  public:
    typedef Segment_3        result_type;

    Rep // Segment_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q) const
    {  return Rep(p, q); }

    Segment_3
    operator()( const Point_3& p, const Point_3& q) const
    { return this->operator()(Return_base_tag(), p, q); }
  };




  template <typename K>
  class Construct_source_2
  {
    typedef typename K::Segment_2  Segment_2;
    typedef typename K::Ray_2      Ray_2;
    typedef typename K::Point_2    Point_2;
  public:
    typedef const Point_2&                result_type;

    result_type
    operator()(const Segment_2& s) const
    {  return s.rep().source(); }

    result_type
    operator()(const Ray_2& r) const
    {  return r.rep().source(); }
  };

  template <typename K>
  class Construct_source_3
  {
    typedef typename K::Segment_3  Segment_3;
    typedef typename K::Ray_3      Ray_3;
    typedef typename K::Point_3    Point_3;
  public:
    typedef const Point_3&         result_type;

    result_type
    operator()(const Segment_3& s) const
    {  return s.rep().source(); }

    result_type
    operator()(const Ray_3& r) const
    {  return r.rep().source(); }
  };


  template <typename K>
  class Construct_target_2
  {
    typedef typename K::Segment_2  Segment_2;
    typedef typename K::Point_2    Point_2;
  public:
    typedef const Point_2&         result_type;

    result_type
    operator()(const Segment_2& s) const
    {  return s.rep().target(); }
  };

  template <typename K>
  class Construct_target_3
  {
    typedef typename K::Segment_3  Segment_3;
    typedef typename K::Point_3    Point_3;
  public:
    typedef const Point_3&         result_type;

    result_type
    operator()(const Segment_3& s) const
    {  return s.rep().target(); }
  };

  template <typename K>
  class Construct_second_point_2
  {
    typedef typename K::Ray_2    Ray_2;
    typedef typename K::Point_2  Point_2;
  public:
    typedef const Point_2&       result_type;

    result_type
    operator()(const Ray_2& r) const
    {  return r.rep().second_point(); }
  };

  template <typename K>
  class Construct_second_point_3
  {
    typedef typename K::Ray_3    Ray_3;
    typedef typename K::Point_3  Point_3;
  public:
    typedef Point_3              result_type;

    result_type // const result_type& // Homogeneous...
    operator()(const Ray_3& r) const
    {  return r.rep().second_point(); }
  };

  template <typename K>
  class Construct_sphere_3
  {
    typedef typename K::FT         FT;
    typedef typename K::Point_3    Point_3;
    typedef typename K::Sphere_3   Sphere_3;
    typedef typename K::Circle_3   Circle_3;
    typedef typename Sphere_3::Rep Rep;
  public:
    typedef Sphere_3               result_type;

    Rep // Sphere_3
    operator()(Return_base_tag, const Point_3& center, const FT& squared_radius,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Rep(center, squared_radius, orientation); }

    Rep // Sphere_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    {  return Rep(p, q, r, s); }

    Rep // Sphere_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q, const Point_3& r,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Rep(p, q, r, orientation); }

    Rep // Sphere_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Rep(p, q, orientation); }

    Rep // Sphere_3
    operator()(Return_base_tag, const Point_3& center,
	        Orientation orientation = COUNTERCLOCKWISE) const
    {  return Rep(center, orientation); }

    Rep
    operator() (Return_base_tag, const Circle_3 & c) const
    { return c.rep().diametral_sphere(); }

    Sphere_3
    operator()( const Point_3& center, const FT& squared_radius,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return this->operator()(Return_base_tag(), center, squared_radius, orientation); }

    Sphere_3
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { return this->operator()(Return_base_tag(), p, q, r, s); }

    Sphere_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return this->operator()(Return_base_tag(), p, q, r, orientation); }

    Sphere_3
    operator()( const Point_3& p, const Point_3& q,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return this->operator()(Return_base_tag(), p, q, orientation); }

    Sphere_3
    operator()( const Point_3& center,
	        Orientation orientation = COUNTERCLOCKWISE) const
    { return this->operator()(Return_base_tag(), center, orientation); }

    Sphere_3
    operator() (const Circle_3 & c) const
    { return this->operator()(Return_base_tag(), c); }

  };

  template <typename K>
  class Construct_supporting_plane_3
  {
    typedef typename K::Triangle_3  Triangle_3;
    typedef typename K::Plane_3     Plane_3;
  public:
    typedef Plane_3          result_type;

    Plane_3
    operator()( const Triangle_3& t) const
    { return t.rep().supporting_plane(); }

  };

  template <typename K>
  class Construct_tetrahedron_3
  {
    typedef typename K::Tetrahedron_3   Tetrahedron_3;
    typedef typename K::Point_3         Point_3;
    typedef typename Tetrahedron_3::Rep Rep;
  public:
    typedef Tetrahedron_3    result_type;

    Rep // Tetrahedron_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { return Rep(p, q, r, s); }

    Tetrahedron_3
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { return this->operator()(Return_base_tag(), p, q, r, s); }
  };

  template <typename K>
  class Construct_triangle_2
  {
    typedef typename K::Triangle_2   Triangle_2;
    typedef typename Triangle_2::Rep  Rep;
    typedef typename K::Point_2      Point_2;
  public:
    typedef Triangle_2       result_type;

    Rep // Triangle_2
    operator()(Return_base_tag, const Point_2& p, const Point_2& q, const Point_2& r) const
    { return Rep(p, q, r); }

    Triangle_2
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    { return this->operator()(Return_base_tag(), p, q, r); }
  };

  template <typename K>
  class Construct_triangle_3
  {
    typedef typename K::Triangle_3   Triangle_3;
    typedef typename K::Point_3      Point_3;
    typedef typename Triangle_3::Rep Rep;
  public:
    typedef Triangle_3       result_type;

    Rep // Triangle_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q, const Point_3& r) const
    { return Rep(p, q, r); }

    Triangle_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& r) const
    { return this->operator()(Return_base_tag(), p, q, r); }
  };

  template <typename K>
  class Construct_unit_normal_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Vector_3         Vector_3;
  public:
    typedef Vector_3           result_type;

    Vector_3
    operator()(const Point_3& p,const Point_3& q, const Point_3& r) const
    { 
      CGAL_kernel_precondition(! K().collinear_3_object()(p,q,r) ); 
      Vector_3 res = CGAL::cross_product(q-p, r-p);
      res = res / CGAL::sqrt(res.squared_length());
      return res; 
	}
  };

  template <typename K>
  class Construct_vertex_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Triangle_3       Triangle_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
  public:
    template<typename>
    struct result {
      typedef const Point_3& type;
    };

    template<typename T>
    struct result<T(Iso_cuboid_3, int)> {
      typedef Point_3 type;
    };

    const Point_3&
    operator()( const Segment_3& s, int i) const
    { return s.rep().vertex(i); }

    const Point_3&
    operator()( const Triangle_3& t, int i) const
    { return t.rep().vertex(i); }

    Point_3
    operator()( const Iso_cuboid_3& r, int i) const
      { return r.rep().vertex(i); }

    const Point_3&
    operator()( const Tetrahedron_3& t, int i) const
    { return t.rep().vertex(i); }
  };

  template <typename K>
  class Construct_cartesian_const_iterator_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Vector_2         Vector_2;
    typedef typename K::Cartesian_const_iterator_2
    Cartesian_const_iterator_2;

  public:
    typedef Cartesian_const_iterator_2 result_type;

    Cartesian_const_iterator_2
    operator()( const Point_2& p) const
    {
      return p.rep().cartesian_begin();
    }

    Cartesian_const_iterator_2
    operator()( const Point_2& p, int) const
    {
      return p.rep().cartesian_end();
    }

    Cartesian_const_iterator_2
    operator()( const Vector_2& v) const
    {
      return v.rep().cartesian_begin();
    }

    Cartesian_const_iterator_2
    operator()( const Vector_2& v, int) const
    {
      return v.rep().cartesian_end();
    }
  };

  template <typename K>
  class Construct_cartesian_const_iterator_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Vector_3         Vector_3;
    typedef typename K::Cartesian_const_iterator_3
    Cartesian_const_iterator_3;

  public:
    typedef Cartesian_const_iterator_3 result_type;

    Cartesian_const_iterator_3
    operator()( const Point_3& p) const
    {
      return p.rep().cartesian_begin();
    }

    Cartesian_const_iterator_3
    operator()( const Point_3& p, int) const
    {
      return p.rep().cartesian_end();
    }

    Cartesian_const_iterator_3
    operator()( const Vector_3& v) const
    {
      return v.rep().cartesian_begin();
    }

    Cartesian_const_iterator_3
    operator()( const Vector_3& v, int) const
    {
      return v.rep().cartesian_end();
    }
  };

  template <typename K>
  class Coplanar_3
  {
    typedef typename K::Point_3       Point_3;
    typedef typename K::Orientation_3 Orientation_3;
    Orientation_3 o;
  public:
    typedef typename K::Boolean       result_type;

    Coplanar_3() {}
    Coplanar_3(const Orientation_3& o_) : o(o_) {}

    result_type
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    {
      return o(p, q, r, s) == COPLANAR;
    }
  };

  template <typename K>
  class Counterclockwise_in_between_2
  {
    typedef typename K::Direction_2  Direction_2;
  public:
    typedef typename K::Boolean      result_type;

    result_type
    operator()( const Direction_2& p, const Direction_2& q,
	        const Direction_2& r) const
    {
        if ( q < p)
            return ( p < r )||( r <= q );
        else
            return ( p < r )&&( r <= q );
    }
  };

  template <typename K>
  class Do_intersect_2
  {
  public:
    typedef typename K::Boolean     result_type;

    // There are 36 combinaisons, so I use a template.
    template <class T1, class T2>
    result_type
    operator()(const T1& t1, const T2& t2) const
    { return internal::do_intersect(t1, t2, K()); }
  };

  template <typename K>
  class Do_intersect_3
  {
  public:
    typedef typename K::Boolean     result_type;

    // There are x combinaisons, so I use a template.
    template <class T1, class T2>
    result_type
    operator()(const T1& t1, const T2& t2) const
    { return internal::do_intersect(t1, t2, K()); }
  };

  template <typename K>
  class Equal_2
  {
    typedef typename K::Point_2       Point_2;
    typedef typename K::Vector_2      Vector_2;
    typedef typename K::Direction_2   Direction_2;
    typedef typename K::Segment_2     Segment_2;
    typedef typename K::Ray_2         Ray_2;
    typedef typename K::Line_2        Line_2;
    typedef typename K::Triangle_2    Triangle_2;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;
    typedef typename K::Circle_2      Circle_2;

  public:
    typedef typename K::Boolean       result_type;

    result_type
    operator()(const Point_2 &p, const Point_2 &q) const
    {
      return p.rep() == q.rep();
    }

    result_type
    operator()(const Vector_2 &v1, const Vector_2 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Vector_2 &v, const Null_vector &n) const
    {
      return v.rep() == n;
    }

    result_type
    operator()(const Direction_2 &d1, const Direction_2 &d2) const
    {
      return d1.rep() == d2.rep();
    }

    result_type
    operator()(const Segment_2 &s1, const Segment_2 &s2) const
    {
      return s1.source() == s2.source() && s1.target() == s2.target();
    }

    result_type
    operator()(const Line_2 &l1, const Line_2 &l2) const
    {
      return l1.rep() == l2.rep();
    }

    result_type
    operator()(const Ray_2& r1, const Ray_2& r2) const
    {
      return r1.source() == r2.source() && r1.direction() == r2.direction();
    }

    result_type
    operator()(const Circle_2& c1, const Circle_2& c2) const
    {
      return c1.center() == c2.center() &&
	c1.squared_radius() == c2.squared_radius() &&
	c1.orientation() == c2.orientation();
    }

    result_type
    operator()(const Triangle_2& t1, const Triangle_2& t2) const
    {
      int i;
      for(i=0; i<3; i++)
	if ( t1.vertex(0) == t2.vertex(i) )
	  break;

      return (i<3) && t1.vertex(1) == t2.vertex(i+1)
                   && t1.vertex(2) == t2.vertex(i+2);
    }

    result_type
    operator()(const Iso_rectangle_2& i1, const Iso_rectangle_2& i2) const
    {
      return ((i1.min)() == (i2.min)()) && ((i1.max)() == (i2.max)());
    }
  };

  template <typename K>
  class Equal_3
  {
    typedef typename K::Point_3       Point_3;
    typedef typename K::Vector_3      Vector_3;
    typedef typename K::Direction_3   Direction_3;
    typedef typename K::Segment_3     Segment_3;
    typedef typename K::Line_3        Line_3;
    typedef typename K::Ray_3         Ray_3;
    typedef typename K::Triangle_3    Triangle_3;
    typedef typename K::Tetrahedron_3 Tetrahedron_3;
    typedef typename K::Sphere_3      Sphere_3;
    typedef typename K::Iso_cuboid_3  Iso_cuboid_3;
    typedef typename K::Plane_3       Plane_3;
    typedef typename K::Circle_3      Circle_3;

  public:
    typedef typename K::Boolean       result_type;

    // Point_3 is special case since the global operator== would recurse.
    result_type
    operator()(const Point_3 &p, const Point_3 &q) const
    {
      return p.x() == q.x() && p.y() == q.y() && p.z() == q.z();
    }

    result_type
    operator()(const Plane_3 &v1, const Plane_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Iso_cuboid_3 &v1, const Iso_cuboid_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Sphere_3 &v1, const Sphere_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Tetrahedron_3 &v1, const Tetrahedron_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Triangle_3 &v1, const Triangle_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Ray_3 &v1, const Ray_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Line_3 &v1, const Line_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Direction_3 &v1, const Direction_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Segment_3 &v1, const Segment_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Vector_3 &v1, const Vector_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }

    result_type
    operator()(const Vector_3 &v, const Null_vector &n) const
    {
      return v.rep() == n;
    }

    result_type
    operator()(const Circle_3 &v1, const Circle_3 &v2) const
    {
      return v1.rep() == v2.rep();
    }
  };

  template <typename K>
  class Has_on_boundary_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_boundary(p); }

    result_type
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_boundary(p); }

    result_type
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return K().bounded_side_2_object()(r,p) == ON_BOUNDARY; }
  };

  template <typename K>
  class Has_on_boundary_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
    typedef typename K::Plane_3          Plane_3;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.rep().has_on_boundary(p); }

    result_type
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.rep().has_on_boundary(p); }

    result_type
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return c.rep().has_on_boundary(p); }

  };

  template <typename K>
  class Has_on_bounded_side_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_bounded_side(p); }

    result_type
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_bounded_side(p); }

    result_type
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { return K().bounded_side_2_object()(r,p) == ON_BOUNDED_SIDE; }
  };

  template <typename K>
  class Has_on_bounded_side_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
    typedef typename K::Circle_3         Circle_3;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.has_on_bounded_side(p); }

    result_type
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.rep().has_on_bounded_side(p); }

    result_type
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return c.rep().has_on_bounded_side(p); }

    result_type
    operator()(const Circle_3& c, const Point_3& p) const
    {
      CGAL_kernel_precondition(
        K().has_on_3_object()(c.supporting_plane(),p)
      );
      return c.rep().has_on_bounded_side(p); 
    }
  };

  template <typename K>
  class Has_on_negative_side_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Line_2           Line_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_negative_side(p); }

    result_type
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_negative_side(p); }

    result_type
    operator()( const Line_2& l, const Point_2& p) const
    { return l.has_on_negative_side(p); }
  };

  template <typename K>
  class Has_on_negative_side_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Plane_3          Plane_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.has_on_negative_side(p); }

    result_type
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.rep().has_on_negative_side(p); }

    result_type
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.rep().has_on_negative_side(p); }
  };

  template <typename K>
  class Has_on_positive_side_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Line_2           Line_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_positive_side(p); }

    result_type
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_positive_side(p); }

    result_type
    operator()( const Line_2& l, const Point_2& p) const
    { return l.has_on_positive_side(p); }
  };

  template <typename K>
  class Has_on_positive_side_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Plane_3          Plane_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.has_on_positive_side(p); }

    result_type
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.rep().has_on_positive_side(p); }

    result_type
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.rep().has_on_positive_side(p); }
  };

  template <typename K>
  class Has_on_unbounded_side_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Circle_2         Circle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Circle_2& c, const Point_2& p) const
    { return c.has_on_unbounded_side(p); }

    result_type
    operator()( const Triangle_2& t, const Point_2& p) const
    { return t.has_on_unbounded_side(p); }

    result_type
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    {
      return K().bounded_side_2_object()(r,p)== ON_UNBOUNDED_SIDE;
    }

  };

  template <typename K>
  class Has_on_unbounded_side_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Circle_3         Circle_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.has_on_unbounded_side(p); }

    result_type
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.rep().has_on_unbounded_side(p); }

    result_type
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return c.rep().has_on_unbounded_side(p); }

    result_type
    operator()(const Circle_3& c, const Point_3& p) const
    {
      CGAL_kernel_precondition(
        K().has_on_3_object()(c.supporting_plane(),p)
      );
      return c.rep().has_on_unbounded_side(p); 
    }
  };

  template <typename K>
  class Has_on_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Line_2           Line_2;
    typedef typename K::Ray_2            Ray_2;
    typedef typename K::Segment_2        Segment_2;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Line_2& l, const Point_2& p) const
    { return l.has_on(p); }

    result_type
    operator()( const Ray_2& r, const Point_2& p) const
    { return r.has_on(p); }

    result_type
    operator()( const Segment_2& s, const Point_2& p) const
    { return s.has_on(p); }
  };

  template <typename K>
  class Intersect_2
  {
    typedef typename K::Object_2    Object_2;
  public:
    typedef Object_2                result_type;

    // 25 possibilities, so I keep the template.
    template <class T1, class T2>
    Object_2
    operator()(const T1& t1, const T2& t2) const
    { return internal::intersection(t1, t2, K()); }
  };

  template <typename K>
  class Intersect_3
  {
    typedef typename K::Object_3    Object_3;
    typedef typename K::Plane_3     Plane_3;
  public:
    typedef Object_3                result_type;

    // n possibilities, so I keep the template.
    template <class T1, class T2>
    Object_3
    operator()(const T1& t1, const T2& t2) const
    { return internal::intersection(t1, t2, K() ); }

    Object_3
    operator()(const Plane_3& pl1, const Plane_3& pl2, const Plane_3& pl3)const
    { return internal::intersection(pl1, pl2, pl3, K() ); }
  };

  template <typename K>
  class Is_degenerate_2
  {
    typedef typename K::Circle_2          Circle_2;
    typedef typename K::Iso_rectangle_2   Iso_rectangle_2;
    typedef typename K::Line_2            Line_2;
    typedef typename K::Ray_2             Ray_2;
    typedef typename K::Segment_2         Segment_2;
    typedef typename K::Triangle_2        Triangle_2;
    typedef typename K::Circle_3          Circle_3;
  public:
    typedef typename K::Boolean           result_type;

    result_type
    operator()( const Circle_2& c) const
    { return c.is_degenerate(); }

    result_type
    operator()( const Iso_rectangle_2& r) const
    { return (r.xmin() == r.xmax()) || (r.ymin() == r.ymax()); }

    result_type
    operator()( const Line_2& l) const
    { return CGAL_NTS is_zero(l.a())  && CGAL_NTS is_zero(l.b()); }

    result_type
    operator()( const Ray_2& r) const
    { return r.rep().is_degenerate(); }

    result_type
    operator()( const Segment_2& s) const
    { return s.source() == s.target(); }

    result_type
    operator()( const Triangle_2& t) const
    { return t.is_degenerate(); }

    result_type
    operator()( const Circle_3& c) const
    { return c.rep().is_degenerate(); }
  };

  template <typename K>
  class Is_degenerate_3
  {
    typedef typename K::Iso_cuboid_3      Iso_cuboid_3;
    typedef typename K::Line_3            Line_3;
    typedef typename K::Circle_3          Circle_3;
    typedef typename K::Plane_3           Plane_3;
    typedef typename K::Ray_3             Ray_3;
    typedef typename K::Segment_3         Segment_3;
    typedef typename K::Sphere_3          Sphere_3;
    typedef typename K::Triangle_3        Triangle_3;
    typedef typename K::Tetrahedron_3     Tetrahedron_3;
  public:
    typedef typename K::Boolean           result_type;

    result_type
    operator()( const Iso_cuboid_3& c) const
    { return c.rep().is_degenerate(); }

    result_type
    operator()( const Line_3& l) const
    { return l.rep().is_degenerate();  }

    result_type
    operator()( const Plane_3& pl) const
    { return pl.rep().is_degenerate(); }

    result_type
    operator()( const Ray_3& r) const
    { return r.rep().is_degenerate(); }

    result_type
    operator()( const Segment_3& s) const
    { return s.rep().is_degenerate(); }

    result_type
    operator()( const Sphere_3& s) const
    { return s.rep().is_degenerate(); }

    result_type
    operator()( const Triangle_3& t) const
    { return t.rep().is_degenerate(); }

    result_type
    operator()( const Tetrahedron_3& t) const
    { return t.rep().is_degenerate(); }

    result_type
    operator()( const Circle_3& t) const
    { return t.rep().is_degenerate(); }

  };

  template <typename K>
  class Is_horizontal_2
  {
    typedef typename K::Line_2    Line_2;
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Ray_2     Ray_2;
  public:
    typedef typename K::Boolean   result_type;

    result_type
    operator()( const Line_2& l) const
    { return CGAL_NTS is_zero(l.a()); }

    result_type
    operator()( const Segment_2& s) const
    { return s.is_horizontal(); }

    result_type
    operator()( const Ray_2& r) const
    { return r.is_horizontal(); }
  };

  template <typename K>
  class Is_vertical_2
  {
    typedef typename K::Line_2    Line_2;
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Ray_2     Ray_2;
  public:
    typedef typename K::Boolean   result_type;

    result_type
    operator()( const Line_2& l) const
    { return CGAL_NTS is_zero(l.b()); }

    result_type
    operator()( const Segment_2& s) const
    { return s.is_vertical(); }

    result_type
    operator()( const Ray_2& r) const
    { return r.is_vertical(); }
  };

  template <typename K>
  class Left_turn_2
  {
    typedef typename K::Point_2        Point_2;
    typedef typename K::Orientation_2  Orientation_2;
    Orientation_2 o;
  public:
    typedef typename K::Boolean        result_type;

    Left_turn_2() {}
    Left_turn_2(const Orientation_2& o_) : o(o_) {}

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return o(p, q, r) == LEFT_TURN; }
  };

  template <typename K>
  class Less_rotate_ccw_2
  {
    typedef typename K::Point_2        Point_2;
    typedef typename K::Orientation_2  Orientation_2;
    typedef typename K::Collinear_are_ordered_along_line_2
    Collinear_are_ordered_along_line_2;
    Orientation_2 o;
    Collinear_are_ordered_along_line_2 co;
  public:
    typedef typename K::Boolean        result_type;

    Less_rotate_ccw_2() {}
    Less_rotate_ccw_2(const Orientation_2& o_,
		      const Collinear_are_ordered_along_line_2& co_)
      : o(o_), co(co_)
    {}

    result_type
    operator()(const Point_2& r, const Point_2& p, const Point_2& q) const
    {
      typename K::Orientation ori = o(r, p, q);
      if ( ori == LEFT_TURN )
	return true;
      else if ( ori == RIGHT_TURN )
	return false;
      else
	{
	  if (p == r) return false;
	  if (q == r) return true;
	  if (p == q) return false;
	  return co( r, q, p);
	}
    }
  };

  template <typename K>
  class Oriented_side_3
  {
    typedef typename K::Point_3        Point_3;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Plane_3        Plane_3;
    typedef typename K::Sphere_3       Sphere_3;
  public:
    typedef typename K::Oriented_side  result_type;

    result_type
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.rep().oriented_side(p); }

    result_type
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.rep().oriented_side(p); }

    result_type
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    { return t.rep().oriented_side(p); }
  };

} // namespace CommonKernelFunctors
} //namespace CGAL

#endif // CGAL_KERNEL_FUNCTION_OBJECTS_H

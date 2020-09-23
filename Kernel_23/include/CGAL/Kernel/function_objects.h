// Copyright (c) 1999,2002,2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany)
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra, Sylvain Pion,
//                 Camille Wormser, Stephane Tayeb, Pierre Alliez



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


#include <cmath> // for Compute_dihedral_angle

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
          // the two cosine are >= 0, square(cosine) is decreasing on [0,pi/2]
          return CGAL::compare(CGAL::square(cosine)*
                               abac1.squared_length()*abad1.squared_length(),
                               CGAL::square(sc_prod_1));
        }
        else {
          return SMALLER;
        }
      }
      else {
        if(cosine < 0) {
          // the two cosine are < 0, square(cosine) is increasing on [pi/2,pi]
          return CGAL::compare(CGAL::square(sc_prod_1),
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
          return CGAL::compare(CGAL::square(sc_prod_2)*
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
          return CGAL::compare(CGAL::square(sc_prod_1)*
                               abac2.squared_length()*abad2.squared_length(),
                               CGAL::square(sc_prod_2)*
                               abac1.squared_length()*abad1.squared_length());
        }
        else
          return LARGER;
        }
    }
  };

  template < typename K >
  class Compare_power_distance_3
  {
  public:
    typedef typename K::Weighted_point_3                  Weighted_point_3;
    typedef typename K::Point_3                           Point_3;
    typedef typename K::Comparison_result                 Comparison_result;

    typedef Comparison_result                             result_type;

    Comparison_result operator()(const Point_3 & p,
                                 const Weighted_point_3 & q,
                                 const Weighted_point_3 & r) const
    {
      return compare_power_distanceC3(p.x(), p.y(), p.z(),
                                      q.x(), q.y(), q.z(), q.weight(),
                                      r.x(), r.y(), r.z(), r.weight());
    }
  };

  template < typename K >
  class Construct_weighted_circumcenter_3
  {
  public:
    typedef typename K::Weighted_point_3               Weighted_point_3;
    typedef typename K::Point_3                        Point_3;
    typedef typename K::FT                             FT;

    typedef Point_3                                    result_type;

    Point_3 operator()(const Weighted_point_3 & p,
                       const Weighted_point_3 & q,
                       const Weighted_point_3 & r,
                       const Weighted_point_3 & s) const
    {
      FT x, y, z;
      weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                              q.x(), q.y(), q.z(), q.weight(),
                              r.x(), r.y(), r.z(), r.weight(),
                              s.x(), s.y(), s.z(), s.weight(),
                              x,y,z);
      return Point_3(x,y,z);
    }

    Point_3 operator()(const Weighted_point_3 & p,
                       const Weighted_point_3 & q,
                       const Weighted_point_3 & r) const
    {
      FT x, y, z;
      weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                              q.x(), q.y(), q.z(), q.weight(),
                              r.x(), r.y(), r.z(), r.weight(),
                              x,y,z);
      return Point_3(x,y,z);
    }

    Point_3 operator()(const Weighted_point_3 & p,
                       const Weighted_point_3 & q) const
    {
      FT x, y, z;
      weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                              q.x(), q.y(), q.z(), q.weight(),
                              x,y,z);
      return Point_3(x,y,z);
    }
  };

  template < class K >
  class Power_side_of_bounded_power_circle_2
  {
  public:
    typedef typename K::Weighted_point_2               Weighted_point_2;
    typedef Bounded_side                               result_type;

    Bounded_side operator()(const Weighted_point_2& p,
                            const Weighted_point_2& q,
                            const Weighted_point_2& r,
                            const Weighted_point_2& t) const
    {
      K traits;
      typename K::Orientation_2 orientation = traits.orientation_2_object();
      typename K::Construct_point_2 wp2p = traits.construct_point_2_object();
      typename K::Power_side_of_oriented_power_circle_2 power_test =
        traits.power_side_of_oriented_power_circle_2_object();
      typename K::Orientation o = orientation(wp2p(p),wp2p(q),wp2p(r));
      typename K::Oriented_side os = power_test(p,q,r,t);

      CGAL_assertion(o != COPLANAR);
      return enum_cast<Bounded_side>(o * os);
    }

    Bounded_side operator()(const Weighted_point_2& p,
                            const Weighted_point_2& q,
                            const Weighted_point_2& t) const
    {
      return power_side_of_bounded_power_circleC2(p.x(), p.y(), p.weight(),
                                                  q.x(), q.y(), q.weight(),
                                                  t.x(), t.y(), t.weight());
    }

    Bounded_side operator()(const Weighted_point_2& p,
                            const Weighted_point_2& t) const
    {
      return enum_cast<Bounded_side>(
            - CGAL_NTS sign( CGAL_NTS square(p.x() - t.x()) +
                             CGAL_NTS square(p.y() - t.y()) +
                             p.weight() - t.weight()) );
    }
  };

  // operator ()
  // return the sign of the power test of  last weighted point
  // with respect to the smallest sphere orthogonal to the others
  template< typename K >
  class Power_side_of_bounded_power_sphere_3
  {
  public:
    typedef typename K::Weighted_point_3               Weighted_point_3;
    typedef typename K::Sign                           Sign;

    typedef Bounded_side                               result_type;

    Bounded_side operator()(const Weighted_point_3 & p,
                            const Weighted_point_3 & q,
                            const Weighted_point_3 & r,
                            const Weighted_point_3 & s,
                            const Weighted_point_3 & t) const
    {
      K traits;
      typename K::Orientation_3  orientation = traits.orientation_3_object();
      typename K::Construct_point_3 wp2p = traits.construct_point_3_object();
      typename K::Power_side_of_oriented_power_sphere_3 power_test =
          traits.power_side_of_oriented_power_sphere_3_object();
      typename K::Orientation o = orientation(wp2p(p),wp2p(q),wp2p(r),wp2p(s));
      typename K::Oriented_side os = power_test(p,q,r,s,t);
      // Power_side_of_oriented_power_sphere_3
      // returns in fact minus the 5x5 determinant of lifted (p,q,r,s,t)
      CGAL_assertion(o != COPLANAR);
      return enum_cast<Bounded_side>(o * os);
    }

    Bounded_side operator()(const Weighted_point_3 & p,
                            const Weighted_point_3 & q,
                            const Weighted_point_3 & r,
                            const Weighted_point_3 & s) const
    {
      return power_side_of_bounded_power_sphereC3(
            p.x(), p.y(), p.z(), p.weight(),
            q.x(), q.y(), q.z(), q.weight(),
            r.x(), r.y(), r.z(), r.weight(),
            s.x(), s.y(), s.z(), s.weight());
    }

    Bounded_side operator()(const Weighted_point_3 & p,
                            const Weighted_point_3 & q,
                            const Weighted_point_3 & r) const
    {
      return power_side_of_bounded_power_sphereC3(
            p.x(), p.y(), p.z(), p.weight(),
            q.x(), q.y(), q.z(), q.weight(),
            r.x(), r.y(), r.z(), r.weight());
    }

    Bounded_side operator()(const Weighted_point_3 & p,
                            const Weighted_point_3 & q) const
    {
      return enum_cast<Bounded_side>(
            - CGAL_NTS sign( CGAL_NTS square(p.x()-q.x()) +
                             CGAL_NTS square(p.y()-q.y()) +
                             CGAL_NTS square(p.z()-q.z()) +
                             p.weight() - q.weight()));
    }
  };

  template < typename K >
  class Power_side_of_oriented_power_sphere_3
  {
  public:
    typedef typename K::Weighted_point_3                  Weighted_point_3;
    typedef typename K::Oriented_side                     Oriented_side;

    typedef Oriented_side                                 result_type;

    Oriented_side operator()(const Weighted_point_3 & p,
                             const Weighted_point_3 & q,
                             const Weighted_point_3 & r,
                             const Weighted_point_3 & s,
                             const Weighted_point_3 & t) const
    {
      return power_side_of_oriented_power_sphereC3(p.x(), p.y(), p.z(), p.weight(),
                                                   q.x(), q.y(), q.z(), q.weight(),
                                                   r.x(), r.y(), r.z(), r.weight(),
                                                   s.x(), s.y(), s.z(), s.weight(),
                                                   t.x(), t.y(), t.z(), t.weight());
    }

    // The methods below are currently undocumented because the definition of
    // orientation is unclear for 3, 2, and 1 point configurations in a 3D space.

    // One should be (very) careful with the order of vertices when using them,
    // as swapping points will change the result and one must therefore have a
    // precise idea of what is the positive orientation in the full space.
    // For example, these functions are (currently) used safely in the regular
    // triangulations classes because we always call them on vertices of
    // triangulation cells, which are always positively oriented.

    Oriented_side operator()(const Weighted_point_3 & p,
                             const Weighted_point_3 & q,
                             const Weighted_point_3 & r,
                             const Weighted_point_3 & s) const
    {
      //CGAL_kernel_precondition( coplanar(p, q, r, s) );
      //CGAL_kernel_precondition( !collinear(p, q, r) );
      return power_side_of_oriented_power_sphereC3(p.x(), p.y(), p.z(), p.weight(),
                                                   q.x(), q.y(), q.z(), q.weight(),
                                                   r.x(), r.y(), r.z(), r.weight(),
                                                   s.x(), s.y(), s.z(), s.weight());
    }

    Oriented_side operator()(const Weighted_point_3 & p,
                             const Weighted_point_3 & q,
                             const Weighted_point_3 & r) const
    {
      //CGAL_kernel_precondition( collinear(p, q, r) );
      //CGAL_kernel_precondition( p.point() != q.point() );
      return power_side_of_oriented_power_sphereC3(p.x(), p.y(), p.z(), p.weight(),
                                                   q.x(), q.y(), q.z(), q.weight(),
                                                   r.x(), r.y(), r.z(), r.weight());
    }

    Oriented_side operator()(const Weighted_point_3 & p,
                             const Weighted_point_3 & q) const
    {
      //CGAL_kernel_precondition( p.point() == r.point() );
      return power_side_of_oriented_power_sphereC3(p.weight(),q.weight());
    }
  };

  template < typename K >
  class Compute_weight_2
  {
  public:
    typedef typename K::Weighted_point_2               Weighted_point_2;
    typedef typename K::FT                             Weight;

    typedef const Weight&     result_type;

    const Weight& operator()(const Weighted_point_2 & p) const
    {
      return p.rep().weight();
    }
  };

  template < typename K >
  class Compute_weight_3
  {
  public:
    typedef typename K::Weighted_point_3               Weighted_point_3;
    typedef typename K::FT                             Weight;

    typedef const Weight&                              result_type;

    const Weight& operator()(const Weighted_point_3 & p) const
    {
      return p.rep().weight();
    }
  };

  template < typename K >
  class Compute_power_product_2
  {
  public:
    typedef typename K::Weighted_point_2               Weighted_point_2;
    typedef typename K::FT                             FT;

    typedef FT                                         result_type;

    FT operator()(const Weighted_point_2 & p,
                  const Weighted_point_2 & q) const
    {
      return power_productC2(p.x(), p.y(), p.weight(),
                             q.x(), q.y(), q.weight());
    }
  };

  template < typename K >
  class Compute_power_product_3
  {
  public:
    typedef typename K::Weighted_point_3               Weighted_point_3;
    typedef typename K::FT                             FT;

    typedef FT                                         result_type;

    FT operator()(const Weighted_point_3 & p,
                  const Weighted_point_3 & q) const
    {
      return power_productC3(p.x(), p.y(), p.z(), p.weight(),
                             q.x(), q.y(), q.z(), q.weight());
    }
  };

  template < typename K >
  class Compute_squared_radius_smallest_orthogonal_circle_2
  {
  public:
    typedef typename K::Weighted_point_2               Weighted_point_2;
    typedef typename K::FT                             FT;

    typedef FT                                         result_type;

    FT operator()(const Weighted_point_2& p,
                  const Weighted_point_2& q,
                  const Weighted_point_2& r) const
    {
      return squared_radius_orthogonal_circleC2(p.x(), p.y(), p.weight(),
                                                q.x(), q.y(), q.weight(),
                                                r.x(), r.y(), r.weight());
    }

    FT operator()(const Weighted_point_2& p,
                  const Weighted_point_2& q) const
    {
      return squared_radius_smallest_orthogonal_circleC2(p.x(), p.y(), p.weight(),
                                                         q.x(), q.y(), q.weight());
    }

    FT operator()(const Weighted_point_2& p) const
    {
      return - p.weight();
    }
  };

  template < typename K >
  class Compute_squared_radius_smallest_orthogonal_sphere_3
  {
  public:
    typedef typename K::Weighted_point_3               Weighted_point_3;
    typedef typename K::FT                             FT;

    typedef FT                                         result_type;

    FT operator()(const Weighted_point_3 & p,
                  const Weighted_point_3 & q,
                  const Weighted_point_3 & r,
                  const Weighted_point_3 & s) const
    {
      return squared_radius_orthogonal_sphereC3(p.x(), p.y(), p.z(), p.weight(),
                                                q.x(), q.y(), q.z(), q.weight(),
                                                r.x(), r.y(), r.z(), r.weight(),
                                                s.x(), s.y(), s.z(), s.weight());
    }

    FT operator()(const Weighted_point_3 & p,
                  const Weighted_point_3 & q,
                  const Weighted_point_3 & r) const
    {
      return squared_radius_smallest_orthogonal_sphereC3(p.x(), p.y(), p.z(), p.weight(),
                                                         q.x(), q.y(), q.z(), q.weight(),
                                                         r.x(), r.y(), r.z(), r.weight());
    }

    FT operator()(const Weighted_point_3 & p,
                  const Weighted_point_3 & q) const
    {
      return squared_radius_smallest_orthogonal_sphereC3(p.x(), p.y(), p.z(), p.weight(),
                                                         q.x(), q.y(), q.z(), q.weight());
    }

    FT operator()(const Weighted_point_3 & p) const
    {
      return - p.weight();
    }
  };

  // Compute the square radius of the sphere centered in t
  // and orthogonal to  the sphere orthogonal to p,q,r,s
  template< typename K>
  class Compute_power_distance_to_power_sphere_3
  {
  public:
    typedef typename K::Weighted_point_3                  Weighted_point_3;
    typedef typename K::FT                                FT;

    typedef FT                                            result_type;

    result_type operator()(const Weighted_point_3 & p,
                           const Weighted_point_3 & q,
                           const Weighted_point_3 & r,
                           const Weighted_point_3 & s,
                           const Weighted_point_3 & t) const
    {
      return power_distance_to_power_sphereC3 (p.x(),p.y(),p.z(),FT(p.weight()),
                                               q.x(),q.y(),q.z(),FT(q.weight()),
                                               r.x(),r.y(),r.z(),FT(r.weight()),
                                               s.x(),s.y(),s.z(),FT(s.weight()),
                                               t.x(),t.y(),t.z(),FT(t.weight()));
    }
  };

  template <typename K>
  class Compare_weighted_squared_radius_3
  {
  public:
    typedef typename K::Weighted_point_3                  Weighted_point_3;
    typedef typename K::Comparison_result                 Comparison_result;
    typedef typename K::FT                                FT;

    typedef Comparison_result                             result_type;

    result_type operator()(const Weighted_point_3 & p,
                           const Weighted_point_3 & q,
                           const Weighted_point_3 & r,
                           const Weighted_point_3 & s,
                           const FT& w) const
    {
      return CGAL::compare(squared_radius_orthogonal_sphereC3(
                             p.x(),p.y(),p.z(),p.weight(),
                             q.x(),q.y(),q.z(),q.weight(),
                             r.x(),r.y(),r.z(),r.weight(),
                             s.x(),s.y(),s.z(),s.weight()),
                           w);
    }

    result_type operator()(const Weighted_point_3 & p,
                           const Weighted_point_3 & q,
                           const Weighted_point_3 & r,
                           const FT& w) const
    {
      return CGAL::compare(squared_radius_smallest_orthogonal_sphereC3(
                             p.x(),p.y(),p.z(),p.weight(),
                             q.x(),q.y(),q.z(),q.weight(),
                             r.x(),r.y(),r.z(),r.weight()),
                           w);
    }

    result_type operator()(const Weighted_point_3 & p,
                           const Weighted_point_3 & q,
                           const FT& w) const
    {
      return CGAL::compare(squared_radius_smallest_orthogonal_sphereC3(
                             p.x(),p.y(),p.z(),p.weight(),
                             q.x(),q.y(),q.z(),q.weight()),
                           w);
    }

    result_type operator()(const Weighted_point_3 & p,
                           const FT& w) const
    {
      return CGAL::compare(-p.weight(), w);
    }
  };

  template <typename K>
  class Compare_slope_3
  {
    typedef typename K::FT                 FT;
    typedef typename K::Point_3 Point_3;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type operator()(const Point_3& p, const Point_3& q, const Point_3& r, const Point_3& s) const
    {
      Comparison_result sign_pq = CGAL::compare(q.z(),p.z());
      Comparison_result sign_rs = CGAL::compare(s.z(),r.z());

      if(sign_pq != sign_rs){
        return CGAL::compare(static_cast<int>(sign_pq), static_cast<int>(sign_rs));
      }

      if((sign_pq == EQUAL) && (sign_rs == EQUAL)){
        return EQUAL;
      }

      CGAL_assertion( (sign_pq == sign_rs) && (sign_pq != EQUAL)  );

      Comparison_result res = CGAL::compare(square(p.z() - q.z()) * (square(r.x()-s.x())+square(r.y()-s.y())),
                                            square(r.z() - s.z()) *  (square(p.x()-q.x())+square(p.y()-q.y())));
      return (sign_pq == SMALLER) ? opposite(res) : res;
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
      return CGAL::compare(squared_distance(p, q), d2);
    }

    template <class T1, class T2, class T3, class T4>
    result_type
    operator()(const T1& p, const T2& q, const T3& r, const T4& s) const
    {
      return CGAL::compare(squared_distance(p, q), squared_distance(r, s));
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
      return CGAL::compare(squared_distance(p, q), d2);
    }

    template <class T1, class T2, class T3, class T4>
    result_type
    operator()(const T1& p, const T2& q, const T3& r, const T4& s) const
    {
      return CGAL::compare(squared_distance(p, q), squared_distance(r, s));
    }
  };

 template <typename K>
 class Compute_approximate_angle_3
 {
   typedef typename K::Point_3 Point_3;
   typedef typename K::Vector_3 Vector_3;

 public:
   typedef typename K::FT       result_type;

    result_type
    operator()(const Vector_3& u, const Vector_3& v) const
   {
     K k;
     typename K::Compute_scalar_product_3 scalar_product =
       k.compute_scalar_product_3_object();

     double product = CGAL::sqrt(to_double(scalar_product(u,u)) * to_double(scalar_product(v,v)));

     if(product == 0)
       return 0;

     // cosine
     double dot = to_double(scalar_product(u,v));
     double cosine = dot / product;

     if(cosine > 1.){
       cosine = 1.;
     }
     if(cosine < -1.){
       cosine = -1.;
     }

     return std::acos(cosine) * 180./CGAL_PI;
   }


   result_type
   operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
   {
     K k;
     typename K::Construct_vector_3 vector = k.construct_vector_3_object();

     Vector_3 u = vector(q,p);
     Vector_3 v = vector(q,r);

     return this->operator()(u,v);
   }
 };

 template <typename K>
 class Compute_approximate_dihedral_angle_3
 {
    typedef typename K::Point_3 Point_3;
 public:
   typedef typename K::FT       result_type;

    result_type
    operator()(const Point_3& a, const Point_3& b, const Point_3& c,  const Point_3& d) const
   {
     K k;
     typename K::Construct_vector_3 vector = k.construct_vector_3_object();
     typename K::Construct_cross_product_vector_3 cross_product =
       k.construct_cross_product_vector_3_object();
     typename K::Compute_squared_distance_3 sq_distance =
       k.compute_squared_distance_3_object();
     typename K::Compute_scalar_product_3 scalar_product =
       k.compute_scalar_product_3_object();

     typedef typename K::Vector_3 Vector_3;
     typedef typename K::FT FT;

     const Vector_3 ab = vector(a,b);
     const Vector_3 ac = vector(a,c);
     const Vector_3 ad = vector(a,d);

     const Vector_3 abad = cross_product(ab,ad);
     const double x = CGAL::to_double(scalar_product(cross_product(ab,ac), abad));
     const double l_ab = CGAL::sqrt(CGAL::to_double(sq_distance(a,b)));
     const double y = l_ab * CGAL::to_double(scalar_product(ac,abad));

     return FT(std::atan2(y, x) * 180 / CGAL_PI );
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
  class Compute_L_infinity_distance_2
  {
    typedef typename K::FT              FT;
    typedef typename K::Point_2         Point_2;

  public:
    typedef FT               result_type;

    result_type
    operator()(const Point_2& p,
               const Point_2& q) const
    {
      return (std::max)( CGAL::abs( K().compute_x_2_object()(p) -  K().compute_x_2_object()(q)),
                         CGAL::abs( K().compute_y_2_object()(p) -  K().compute_y_2_object()(q)) );
    }
  };

  template <typename K>
  class Compute_L_infinity_distance_3
  {
    typedef typename K::FT              FT;
    typedef typename K::Point_3         Point_3;

  public:
    typedef FT               result_type;

    result_type
    operator()(const Point_3& p,
               const Point_3& q) const
    {
      return (std::max)( CGAL::abs( K().compute_x_3_object()(p) -  K().compute_x_3_object()(q)),
                         (std::max)(CGAL::abs( K().compute_y_3_object()(p) -  K().compute_y_3_object()(q)),
                                    CGAL::abs( K().compute_z_3_object()(p) -  K().compute_z_3_object()(q))));
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
    operator()(        const Point_3& p1, const Point_3& p2, const Point_3& p3) const
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
  class Construct_line_line_intersection_point_3
  {
    typedef typename K::Line_3 Line;
    typedef typename K::Point_3 Point;
    typename K::Construct_line_3 construct_line;
  public:
    typedef Point result_type;

    Point
    operator()(const Point& l11, const Point& l12,
               const Point& l21, const Point& l22) const
    {
      Line l1 = construct_line(l11, l12);
      Line l2 = construct_line(l21, l22);

      typename cpp11::result_of<typename K::Intersect_3(Line,Line)>::type
        res = typename K::Intersect_3()(l1,l2);
      CGAL_assertion(res!=boost::none);
      const Point* e_pt = boost::get<Point>(&(*res));
      CGAL_assertion(e_pt!=nullptr);
      return *e_pt;
    }
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
  class Construct_plane_line_intersection_point_3
  {
    typedef typename K::Plane_3 Plane;
    typedef typename K::Line_3 Line;
    typedef typename K::Point_3 Point;
    typename K::Construct_plane_3 construct_plane;
    typename K::Construct_line_3 construct_line;
  public:
    typedef Point result_type;

    Point
    operator()(const Point& p1, const Point& p2, const Point& p3,
               const Point& l1, const Point& l2) const
    {
      Plane plane = construct_plane(p1, p2, p3);
      Line line = construct_line( l1, l2 );

      typename cpp11::result_of<typename K::Intersect_3(Plane,Line)>::type
        res = typename K::Intersect_3()(plane,line);
      CGAL_assertion(res!=boost::none);
      const Point* e_pt = boost::get<Point>(&(*res));
      CGAL_assertion(e_pt!=nullptr);
      return *e_pt;
    }

    Point
    operator()(const Plane& plane,
               const Point& l1, const Point& l2) const
    {
      Line line = construct_line( l1, l2 );

      typename cpp11::result_of<typename K::Intersect_3(Plane,Line)>::type
        res = typename K::Intersect_3()(plane,line);
      CGAL_assertion(res!=boost::none);
      const Point* e_pt = boost::get<Point>(&(*res));
      CGAL_assertion(e_pt!=nullptr);
      return *e_pt;
    }
  };

  template <typename K>
  class Construct_point_on_2
  {
    typedef typename K::FT         FT;
    typedef typename K::Point_2    Point_2;
    typedef typename K::Segment_2  Segment_2;
    typedef typename K::Line_2     Line_2;
    typedef typename K::Ray_2      Ray_2;
  public:
    typedef Point_2          result_type;

    Point_2
    operator()( const Line_2& l, const FT i) const
    { return l.point(i); }

    Point_2
    operator()( const Segment_2& s, int i) const
    { return s.point(i); }

    Point_2
    operator()( const Ray_2& r, const FT i) const
    { return r.point(i); }
  };

  template <typename K>
  class Construct_point_on_3
  {
    typedef typename K::FT         FT;
    typedef typename K::Point_3    Point_3;
    typedef typename K::Segment_3  Segment_3;
    typedef typename K::Line_3     Line_3;
    typedef typename K::Ray_3      Ray_3;
    typedef typename K::Plane_3    Plane_3;
  public:
    typedef Point_3          result_type;

    Point_3
    operator()( const Line_3& l, const FT i) const
    { return l.rep().point(i); }

    Point_3
    operator()( const Segment_3& s, int i) const
    { return s.point(i); }

    Point_3
    operator()( const Ray_3& r, const FT i) const
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
  class Construct_projected_point_3
  {
    bool
    is_inside_triangle_3_aux(const typename K::Vector_3& w,
                             const typename K::Point_3& p1,
                             const typename K::Point_3& p2,
                             const typename K::Point_3& q,
                             typename K::Point_3& result,
                             bool& outside,
                             const K& k)
    {
      typedef typename K::Vector_3 Vector_3;
      typedef typename K::FT FT;

      typename K::Construct_vector_3 vector =
        k.construct_vector_3_object();
      typename K::Construct_projected_point_3 projection =
        k.construct_projected_point_3_object();
      typename K::Construct_line_3 line =
        k.construct_line_3_object();
      typename K::Compute_scalar_product_3 scalar_product =
        k.compute_scalar_product_3_object();
      typename K::Construct_cross_product_vector_3 cross_product =
        k.construct_cross_product_vector_3_object();

      const Vector_3 v = cross_product(vector(p1,p2), vector(p1,q));
      if ( scalar_product(v,w) < FT(0))
      {
        if (   scalar_product(vector(p1,q), vector(p1,p2)) >= FT(0)
            && scalar_product(vector(p2,q), vector(p2,p1)) >= FT(0) )
        {
          result = projection(line(p1, p2), q);
          return true;
        }
        outside = true;
      }

      return false;
    }


    /**
     * Returns the nearest point of p1,p2,p3 from origin
     * @param origin the origin point
     * @param p1 the first point
     * @param p2 the second point
     * @param p3 the third point
     * @param k the kernel
     * @return the nearest point from origin
     */
    typename K::Point_3
    nearest_point_3(const typename K::Point_3& origin,
                    const typename K::Point_3& p1,
                    const typename K::Point_3& p2,
                    const typename K::Point_3& p3,
                    const K& k)
    {
      typedef typename K::FT FT;

      typename K::Compute_squared_distance_3 sq_distance =
        k.compute_squared_distance_3_object();

      const FT dist_origin_p1 = sq_distance(origin,p1);
      const FT dist_origin_p2 = sq_distance(origin,p2);
      const FT dist_origin_p3 = sq_distance(origin,p3);

      if (   dist_origin_p2 >= dist_origin_p1
          && dist_origin_p3 >= dist_origin_p1 )
      {
        return p1;
      }
      if ( dist_origin_p3 >= dist_origin_p2 )
      {
        return p2;
      }

      return p3;
    }

    /**
     * @brief returns true if p is inside triangle t. If p is not inside t,
     * result is the nearest point of t from p. WARNING: it is assumed that
     * t and p are on the same plane.
     * @param p the reference point
     * @param t the triangle
     * @param result if p is not inside t, the nearest point of t from p
     * @param k the kernel
     * @return true if p is inside t
     */
    bool
    is_inside_triangle_3(const typename K::Point_3& p,
                         const typename K::Triangle_3& t,
                         typename K::Point_3& result,
                         const K& k)
    {
      typedef typename K::Point_3 Point_3;
      typedef typename K::Vector_3 Vector_3;

      typename K::Construct_vector_3 vector =
        k.construct_vector_3_object();
      typename K::Construct_vertex_3 vertex_on =
        k.construct_vertex_3_object();
      typename K::Construct_cross_product_vector_3 cross_product =
        k.construct_cross_product_vector_3_object();

      const Point_3& t0 = vertex_on(t,0);
      const Point_3& t1 = vertex_on(t,1);
      const Point_3& t2 = vertex_on(t,2);

      Vector_3 w = cross_product(vector(t0,t1), vector(t1,t2));

      bool outside = false;
      if (   is_inside_triangle_3_aux(w, t0, t1, p, result, outside, k)
          || is_inside_triangle_3_aux(w, t1, t2, p, result, outside, k)
          || is_inside_triangle_3_aux(w, t2, t0, p, result, outside, k) )
      {
        return false;
      }

      if ( outside )
      {
        result = nearest_point_3(p,t0,t1,t2,k);
        return false;
      }
      else
      {
        return true;
      }
    }

    /**
    * @brief returns true if p is inside segment s. If p is not inside s,
    * result is the nearest point of s from p. WARNING: it is assumed that
    * t and p are on the same line.
    * @param query the query point
    * @param s the segment
    * @param closest_point_on_segment if query is not inside s, the nearest point of s from p
    * @param k the kernel
    * @return true if p is inside s
    */
    bool
    is_inside_segment_3(const typename K::Point_3& query,
                        const typename K::Segment_3 & s,
                        typename K::Point_3& closest_point_on_segment,
                        const K& k)
    {
      typename K::Construct_vector_3 vector =
        k.construct_vector_3_object();
      typename K::Construct_vertex_3 vertex_on =
        k.construct_vertex_3_object();
      typename K::Compute_scalar_product_3 scalar_product =
        k.compute_scalar_product_3_object();

      typedef typename K::FT FT;
      typedef typename K::Point_3 Point;

      const Point& a = vertex_on(s, 0);
      const Point& b = vertex_on(s, 1);
      if( scalar_product(vector(a,b), vector(a, query)) < FT(0) )
      {
        closest_point_on_segment = a;
        return false;
      }
      if( scalar_product(vector(b,a), vector(b, query)) < FT(0) )
      {
        closest_point_on_segment = b;
        return false;
      }

      // query is on segment
      return true;
    }

  public:
    typename K::Point_3
    operator()(const typename K::Point_3& origin,
               const typename K::Triangle_3& triangle,
               const K& k)
    {
      typedef typename K::Point_3 Point_3;

      typename K::Construct_supporting_plane_3 supporting_plane =
        k.construct_supporting_plane_3_object();
      typename K::Construct_projected_point_3 projection =
        k.construct_projected_point_3_object();
      typename K::Is_degenerate_3 is_degenerate = k.is_degenerate_3_object();

      const typename K::Plane_3 plane = supporting_plane(triangle);
      if(is_degenerate(plane)) {
        // If the plane is degenerate, then the triangle is degenerate, and
        // one tries to find to which segment it is equivalent.
        typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();
        typename K::Construct_vector_3 vector = k.construct_vector_3_object();
        typename K::Compute_x_3 x = k.compute_x_3_object();
        typename K::Compute_y_3 y = k.compute_y_3_object();
        typename K::Compute_z_3 z = k.compute_z_3_object();
        typedef typename K::FT FT;
        typedef typename K::Vector_3 Vector_3;

        const Point_3& a = vertex(triangle, 0);
        const Point_3& b = vertex(triangle, 1);
        const Point_3& c = vertex(triangle, 2);
        const Vector_3 ab = vector(a, b);
        const Vector_3 ac = vector(a, c);
        const Vector_3 bc = vector(b, c);
        const FT linf_ab = (std::max)((std::max)(x(ab), y(ab)), z(ab));
        const FT linf_ac = (std::max)((std::max)(x(ac), y(ac)), z(ac));
        const FT linf_bc = (std::max)((std::max)(x(bc), y(bc)), z(bc));

        typename K::Construct_segment_3 seg = k.construct_segment_3_object();
        if(linf_ab > linf_ac) {
          if(linf_ab > linf_bc) {
            // ab is the maximal segment
            return this->operator()(origin, seg(a, b), k);
          } else {
            // ab > ac, bc >= ab, use bc
            return this->operator()(origin, seg(b, c), k);
          }
        } else { // ab <= ac
          if(linf_ac > linf_bc) {
            // ac is the maximal segment
            return this->operator()(origin, seg(a, c), k);
          } else {
            // ab <= ac, ac <= bc, use bc
            return this->operator()(origin, seg(b, c), k);
          }
        }
      } // degenerate plane

      // Project origin on triangle supporting plane
      const Point_3 proj = projection(plane, origin);


      Point_3 moved_point;
      bool inside = is_inside_triangle_3(proj,triangle,moved_point,k);

      // If proj is inside triangle, return it
      if ( inside )
      {
        return proj;
      }

      // Else return the constructed point
      return moved_point;
    }

    typename K::Point_3
    operator()(const typename K::Point_3& query,
               const typename K::Segment_3& segment,
               const K& k)
    {

      typename K::Is_degenerate_3 is_degenerate =
          k.is_degenerate_3_object();
      typename K::Construct_vertex_3 vertex =
          k.construct_vertex_3_object();

      if(is_degenerate(segment))
        return vertex(segment, 0);

      if(segment.to_vector() * (query-segment.source()) <= 0)
        return segment.source();
      if(segment.to_vector() * (query-segment.target()) >= 0)
        return segment.target();
      // If proj is inside segment, returns it
      return k.construct_projected_point_3_object()(segment.supporting_line(), query);
    }

    typename K::Point_3
    operator()(const typename K::Point_3& query,
               const typename K::Ray_3& ray,
               const K& k)
    {
      if ( ray.to_vector() * (query-ray.source()) <= 0)
        return ray.source();
      else
      {
        return k.construct_projected_point_3_object()(ray.supporting_line(), query);
      }
    }

    // code for operator for plane and point is defined in
    // CGAL/Cartesian/function_objects.h and CGAL/Homogeneous/function_objects.h
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
    { return Intersections::internal::do_intersect(t1, t2, K()); }
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
    { return Intersections::internal::do_intersect(t1, t2, K()); }

    result_type
    operator()(const typename K::Plane_3& pl1, const typename K::Plane_3& pl2, const typename K::Plane_3& pl3) const
    { return Intersections::internal::do_intersect(pl1, pl2, pl3, K() ); }

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

    bool operator()(const Sphere_3& s1, const Sphere_3& s2,
                    const Point_3& a, const Point_3& b) const
    {
      typedef typename K::Circle_3    Circle_3;
      typedef typename K::Point_3     Point_3;
      typedef typename K::Segment_3   Segment_3;
      typedef typename K::Plane_3     Plane_3;
      typedef typename K::Intersect_3 Intersect_3;

      const Has_on_bounded_side_3& has_on_bounded_side = *this;

      const bool a_in_s1 = has_on_bounded_side(s1, a);
      const bool a_in_s2 = has_on_bounded_side(s2, a);

      if(!(a_in_s1 || a_in_s2)) return false;

      const bool b_in_s1 = has_on_bounded_side(s1, b);
      const bool b_in_s2 = has_on_bounded_side(s2, b);

      if(!(b_in_s1 || b_in_s2)) return false;

      if(a_in_s1 && b_in_s1) return true;
      if(a_in_s2 && b_in_s2) return true;

      if(!K().do_intersect_3_object()(s1, s2)) return false;
      const Circle_3 circ(s1, s2);
      const Plane_3& plane = circ.supporting_plane();
      typename CGAL::cpp11::result_of<Intersect_3(Plane_3, Segment_3)>::type
        optional = K().intersect_3_object()(plane, Segment_3(a, b));
      CGAL_kernel_assertion_msg(bool(optional) == true,
                                "the segment does not intersect the supporting"
                                " plane");
      const Point_3* p = boost::get<Point_3>(&*optional);
      CGAL_kernel_assertion_msg(p != 0,
                                "the segment intersection with the plane is "
                                "not a point");
      return squared_distance(circ.center(), *p) < circ.squared_radius();
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
  public:
    template<typename>
    struct result;

    template<typename F, typename A, typename B>
    struct result<F(A,B)> {
      typedef typename Intersection_traits<K, A, B>::result_type type;
    };

    // 25 possibilities, so I keep the template.
    template <class T1, class T2>
    typename Intersection_traits<K, T1, T2>::result_type
    operator()(const T1& t1, const T2& t2) const
    { return Intersections::internal::intersection(t1, t2, K()); }
  };

  template <typename K>
  class Intersect_3
  {
    typedef typename K::Plane_3     Plane_3;
  public:
    template<typename>
    struct result;

    template<typename F, typename A, typename B>
    struct result<F(A, B)> {
      typedef typename Intersection_traits<K, A, B>::result_type type;
    };

    template<typename F>
    struct result<F(Plane_3, Plane_3, Plane_3)> {
      typedef boost::optional<
        boost::variant< typename K::Point_3,
                        typename K::Line_3,
                        typename K::Plane_3 > > type;
    };

    // n possibilities, so I keep the template.
    template <class T1, class T2>
    typename cpp11::result_of< Intersect_3(T1, T2) >::type
    operator()(const T1& t1, const T2& t2) const
    { return Intersections::internal::intersection(t1, t2, K() ); }

    typename boost::optional< boost::variant< typename K::Point_3, typename K::Line_3, typename K::Plane_3 > >
    operator()(const Plane_3& pl1, const Plane_3& pl2, const Plane_3& pl3)const
    { return Intersections::internal::intersection(pl1, pl2, pl3, K() ); }
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


template < typename K >
class Construct_weighted_circumcenter_2
{
public:
  typedef typename K::Weighted_point_2         Weighted_point_2;
  typedef typename K::Point_2                  Point_2;
  typedef typename K::FT                       FT;

  typedef Point_2       result_type;

  result_type operator() (const Weighted_point_2 & p,
                          const Weighted_point_2 & q,
                          const Weighted_point_2 & r) const
  {
    CGAL_kernel_precondition( ! collinear(p.point(), q.point(), r.point()) );
    FT x,y;
    weighted_circumcenterC2(p.x(),p.y(),p.weight(),
                            q.x(),q.y(),q.weight(),
                            r.x(),r.y(),r.weight(),x,y);
    return Point_2(x,y);
  }
};



} // namespace CommonKernelFunctors
} //namespace CGAL

#endif // CGAL_KERNEL_FUNCTION_OBJECTS_H

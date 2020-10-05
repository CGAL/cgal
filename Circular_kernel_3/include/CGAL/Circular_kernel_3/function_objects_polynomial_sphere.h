// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado,
//             Sebastien Loriot

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_SPHERE_H
#define CGAL_SPHERICAL_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_SPHERE_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/kernel_basic.h>
#include <CGAL/is_iterator.h>

#include <CGAL/Spherical_kernel_intersections.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_circular_arc_point_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_sphere_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_line_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_plane_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_circle_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_line_arc_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_circular_arc_3.h>
#include <CGAL/Circular_kernel_3/internal_function_has_on_spherical_kernel.h>
#include <CGAL/Circular_kernel_3/internal_function_compare_spherical_kernel.h>
#include <CGAL/Circular_kernel_3/internal_function_compare_to_right_spherical_kernel.h>
#include <CGAL/Circular_kernel_3/Intersection_traits.h>

namespace CGAL {

namespace SphericalFunctors {

#define CGAL_SPHERICAL_KERNEL_MACRO_FUNCTOR_COMPARE_(V)\
template < class SK > \
  class Compare_ ##V## _3: public SK::Linear_kernel::Compare_ ##V## _3{\
    typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;\
    typedef typename SK::Point_3 Point_3;\
  public:\
    typedef  typename SK::Linear_kernel::Compare_ ##V## _3::result_type result_type;\
    using SK::Linear_kernel::Compare_ ##V## _3::operator();\
    result_type\
    operator() (const Circular_arc_point_3 &p0,\
                const Circular_arc_point_3 &p1) const\
    { return SphericalFunctors::compare_ ##V <SK>(p0, p1); }\
    result_type\
    operator() (const Circular_arc_point_3 &p0,\
                const Point_3 &p1) const\
    { return SphericalFunctors::compare_ ##V <SK>(p0, p1); }\
    result_type\
    operator() (const Point_3 &p0,\
                const Circular_arc_point_3 &p1) const\
    { return SphericalFunctors::compare_ ##V <SK>(p0, p1); }\
  };\

  CGAL_SPHERICAL_KERNEL_MACRO_FUNCTOR_COMPARE_(x)
  CGAL_SPHERICAL_KERNEL_MACRO_FUNCTOR_COMPARE_(y)
  CGAL_SPHERICAL_KERNEL_MACRO_FUNCTOR_COMPARE_(z)
  CGAL_SPHERICAL_KERNEL_MACRO_FUNCTOR_COMPARE_(xy)
  CGAL_SPHERICAL_KERNEL_MACRO_FUNCTOR_COMPARE_(xyz)

  template <class SK>
  class Compute_circular_x_3
  {
    typedef typename SK::Circular_arc_point_3   Circular_arc_point_3;
    typedef typename SK::Root_of_2              Root_of_2;

  public:

    typedef const Root_of_2&                    result_type;

    result_type operator() (const Circular_arc_point_3 & a) const
    { return (a.rep().x()); }
  };

  template <class SK>
  class Compute_circular_y_3
  {
    typedef typename SK::Circular_arc_point_3   Circular_arc_point_3;
    typedef typename SK::Root_of_2                 Root_of_2;

  public:

    typedef  const Root_of_2&  result_type;

    result_type operator() (const Circular_arc_point_3 & a) const
    { return (a.rep().y()); }
  };

  template <class SK>
  class Compute_circular_z_3
  {
    typedef typename SK::Circular_arc_point_3   Circular_arc_point_3;
    typedef typename SK::Root_of_2                 Root_of_2;

  public:

    typedef const Root_of_2&  result_type;

    result_type operator() (const Circular_arc_point_3 & a) const
    { return (a.rep().z()); }
  };

  template < class SK >
  class Equal_3
    : public SK::Linear_kernel::Equal_3
  {
    typedef typename SK::Linear_kernel LK;
    typedef typename LK::Equal_3 LK_Equal_3;

    typedef typename SK::Point_3 Point_3;
    typedef typename SK::Vector_3 Vector_3;
    typedef typename SK::Direction_3 Direction_3;
    typedef typename SK::Line_3 Line_3;
    typedef typename SK::Ray_3 Ray_3;
    typedef typename SK::Segment_3 Segment_3;
    typedef typename SK::Triangle_3 Triangle_3;
    typedef typename SK::Tetrahedron_3 Tetrahedron_3;
    typedef typename SK::Iso_cuboid_3 Iso_cuboid_3;
    typedef typename SK::Plane_3 Plane_3;
    typedef typename SK::Circular_arc_point_3     Circular_arc_point_3;
    typedef typename SK::Circle_3                 Circle_3;
    typedef typename SK::Sphere_3                 Sphere_3;

    typedef typename SK::Line_arc_3               Line_arc_3;
    typedef typename SK::Circular_arc_3           Circular_arc_3;

  public:

    typedef typename SK::Linear_kernel::Equal_3::result_type result_type;

    using SK::Linear_kernel::Equal_3::operator();

    result_type
    operator() (const Circular_arc_point_3 &c0,
                const Circular_arc_point_3 &c1) const
    { return SphericalFunctors::equal<SK>(c0, c1); }

    result_type
    operator() (const Circular_arc_point_3 &c0,
                const Point_3 &c1) const
    { return SphericalFunctors::equal<SK>(c0, Circular_arc_point_3(c1)); }

    result_type
    operator() (const Point_3 &c0,
                const Circular_arc_point_3 &c1) const
    { return SphericalFunctors::equal<SK>(Circular_arc_point_3(c0), c1); }

    // Our Line_arc_3 dont have orientation
    result_type
    operator() (const Line_arc_3 &l0,
                const Line_arc_3 &l1) const
    { return SphericalFunctors::equal<SK>(l0, l1); }

    // Our Circular_arc_3 dont have orientation (as parameter)
    result_type
    operator() (const Circular_arc_3 &c0,
                const Circular_arc_3 &c1) const
    { return SphericalFunctors::equal<SK>(c0, c1); }

  };

  template < class SK >
  class Construct_circular_arc_point_3
  {
    typedef typename SK::Point_3                            Point_3;
    typedef typename SK::Plane_3                            Plane_3;
    typedef typename SK::Line_3                             Line_3;
    typedef typename SK::Circle_3                           Circle_3;
    typedef typename SK::Sphere_3                           Sphere_3;
    typedef typename SK::Circular_arc_point_3               Circular_arc_point_3;
    typedef typename SK::Kernel_base::Circular_arc_point_3  RCircular_arc_point_3;
    typedef typename SK::Root_of_2                          Root_of_2;
    typedef typename Circular_arc_point_3::Rep              Rep;
    typedef typename Circular_arc_point_3::Root_for_spheres_2_3  Root_for_spheres_2_3;

  public:
    typedef  Circular_arc_point_3 result_type;

    result_type
    operator()(void)
    { return Rep(); }

    result_type
      operator()(const Root_of_2 & x,
                 const Root_of_2 & y,
                 const Root_of_2 & z
                 ) const
    { return Rep(x,y,z); }

    result_type
    operator()(const Root_for_spheres_2_3 & np) const
    { return Rep(np); }

    result_type
    operator()(const Point_3 & p) const
    { return Rep(p); }

    // Not Documented
    result_type
    operator()(const Sphere_3 & s1,
               const Sphere_3 & s2,
               const Sphere_3 & s3,
               const bool less_xyz = true) const
    { return Rep(s1,s2,s3,less_xyz); }

    // Not Documented
    result_type
    operator()(const Plane_3 & p,
               const Sphere_3 & s1,
               const Sphere_3 & s2,
               const bool less_xyz = true) const
    { return Rep(p,s1,s2,less_xyz); }

    // Not Documented
    result_type
    operator()(const Sphere_3 & s1,
               const Plane_3 & p,
               const Sphere_3 & s2,
               const bool less_xyz = true) const
    { return Rep(p,s1,s2,less_xyz); }

    // Not Documented
    result_type
    operator()(const Sphere_3 & s1,
               const Sphere_3 & s2,
               const Plane_3 & p,
               const bool less_xyz = true) const
    { return Rep(p,s1,s2,less_xyz); }

    // Not Documented
    result_type
    operator()(const Plane_3 & p1,
               const Plane_3 & p2,
               const Sphere_3 & s,
               const bool less_xyz = true) const
    { return Rep(p1,p2,s,less_xyz); }

    // Not Documented
    result_type
    operator()(const Plane_3 & p1,
               const Sphere_3 & s,
               const Plane_3 & p2,
               const bool less_xyz = true) const
    { return Rep(p1,p2,s,less_xyz); }

    // Not Documented
    result_type
    operator()(const Sphere_3 & s,
               const Plane_3 & p1,
               const Plane_3 & p2,
               const bool less_xyz = true) const
    { return Rep(p1,p2,s,less_xyz); }

    // Not Documented
    result_type
    operator()(const Line_3 & l,
               const Sphere_3 & s,
               const bool less_xyz = true) const
    { return Rep(l,s,less_xyz); }

    // Not Documented
    result_type
    operator()(const Sphere_3 & s,
               const Line_3 & l,
               const bool less_xyz = true) const
    { return Rep(l,s,less_xyz); }

    // Not Documented
    result_type
    operator()(const Circle_3 & c,
               const Sphere_3 & s,
               const bool less_xyz = true) const
    { return Rep(c,s,less_xyz); }

    // Not Documented
    result_type
    operator()(const Sphere_3 & s,
               const Circle_3 & c,
               const bool less_xyz = true) const
    { return Rep(c,s,less_xyz); }

    // Not Documented
    result_type
    operator()(const Circle_3 & c,
               const Plane_3 & p,
               const bool less_xyz = true) const
    { return Rep(c,p,less_xyz); }

    // Not Documented
    result_type
    operator()(const Plane_3 & p,
               const Circle_3 & c,
               const bool less_xyz = true) const
    { return Rep(c,p,less_xyz); }

  };

  template < class SK >
  class Construct_sphere_3
  {
    typedef typename SK::Circular_arc_3 Circular_arc_3;
    typedef typename SK::Linear_kernel LK;
    typedef typename LK::Point_3 Point_3;
    typedef typename LK::Circle_3 Circle_3;
    typedef typename LK::Sphere_3 Sphere_3;
    typedef typename LK::Construct_sphere_3 LK_Construct_sphere_3;
    typedef typename LK::FT FT;

  public:

    typedef typename SK::Linear_kernel::Construct_sphere_3::result_type result_type;

    result_type
    operator()( Return_base_tag tag, const Point_3& center, const FT& squared_radius,
                Orientation orientation = COUNTERCLOCKWISE) const
    { return LK_Construct_sphere_3()(tag, center, squared_radius, orientation); }

    result_type
    operator()( Return_base_tag tag, const Point_3& p, const Point_3& q,
                const Point_3& r, const Point_3& s) const
    { return LK_Construct_sphere_3()(tag, p, q, r, s); }

    result_type
    operator()( Return_base_tag tag, const Point_3& p, const Point_3& q, const Point_3& r,
                Orientation orientation  = COUNTERCLOCKWISE) const
    { return LK_Construct_sphere_3()(tag, p, q, r, orientation); }

    result_type
    operator()( Return_base_tag tag, const Point_3& p, const Point_3& q,
                Orientation orientation = COUNTERCLOCKWISE) const
    { return LK_Construct_sphere_3()(tag, p, q, orientation); }

    result_type
    operator()( Return_base_tag tag, const Point_3& center,
                Orientation orientation = COUNTERCLOCKWISE) const
    { return LK_Construct_sphere_3()(tag, center, orientation); }

    result_type
    operator() (Return_base_tag tag, const Circle_3 & c) const
    { return LK_Construct_sphere_3()(tag, c); }




    result_type
    operator()( const Point_3& center, const FT& squared_radius,
                Orientation orientation = COUNTERCLOCKWISE) const
    { return LK_Construct_sphere_3()(center, squared_radius, orientation); }

    result_type
    operator()( const Point_3& p, const Point_3& q,
                const Point_3& r, const Point_3& s) const
    { return LK_Construct_sphere_3()(p, q, r, s); }

    result_type
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
                Orientation orientation = COUNTERCLOCKWISE) const
    { return LK_Construct_sphere_3()(p, q, r, orientation); }

    result_type
    operator()( const Point_3& p, const Point_3& q,
                Orientation orientation = COUNTERCLOCKWISE) const
    { return LK_Construct_sphere_3()(p, q, orientation); }

    result_type
    operator()( const Point_3& center,
                Orientation orientation = COUNTERCLOCKWISE) const
    { return LK_Construct_sphere_3()(center, orientation); }

    result_type
    operator() (const Circle_3 & c) const
    { return LK_Construct_sphere_3()(c); }

    result_type
    operator() ( const typename SK::Polynomial_for_spheres_2_3 &eq )
    { return SphericalFunctors::construct_sphere_3<SK>(eq); }

    result_type operator() (const Circular_arc_3 & c) const
    { return c.rep().diametral_sphere(); }

  };

  template < class SK >
  class Construct_plane_3
  {
    typedef typename SK::Circular_arc_3 Circular_arc_3;
    typedef typename SK::Linear_kernel LK;
    typedef typename LK::Construct_plane_3 LK_Construct_plane_3;
    typedef typename LK::RT           RT;
    typedef typename LK::Point_3      Point_3;
    typedef typename LK::Vector_3     Vector_3;
    typedef typename LK::Direction_3  Direction_3;
    typedef typename LK::Line_3       Line_3;
    typedef typename LK::Ray_3        Ray_3;
    typedef typename LK::Segment_3    Segment_3;
    typedef typename LK::Plane_3      Plane_3;
    typedef typename LK::Circle_3     Circle_3;

  public:

    typedef typename SK::Linear_kernel::Construct_plane_3::result_type result_type;

  public:

    result_type
    operator()(Return_base_tag tag, const RT& a, const RT& b, const RT& c, const RT& d) const
    { return LK_Construct_plane_3()(tag, a, b, c, d); }

    result_type
    operator()(Return_base_tag tag, const Point_3& p, const Point_3& q, const Point_3& r) const
    { return LK_Construct_plane_3()(tag, p, q, r); }

    result_type
    operator()(Return_base_tag tag, const Point_3& p, const Direction_3& d) const
    { return LK_Construct_plane_3()(tag, p, d); }

    result_type
    operator()(Return_base_tag tag, const Point_3& p, const Vector_3& v) const
    { return LK_Construct_plane_3()(tag, p, v); }

    result_type
    operator()(Return_base_tag tag, const Line_3& l, const Point_3& p) const
    { return LK_Construct_plane_3()(tag, l, p); }

    result_type
    operator()(Return_base_tag tag, const Ray_3& r, const Point_3& p) const
    { return LK_Construct_plane_3()(tag, r, p); }

    result_type
    operator()(Return_base_tag tag, const Segment_3& s, const Point_3& p) const
    { return LK_Construct_plane_3()(tag, s, p); }

    result_type
    operator()(Return_base_tag tag, const Circle_3 & c) const
    { return LK_Construct_plane_3()(tag, c); }

    result_type
    operator()(const RT& a, const RT& b, const RT& c, const RT& d) const
    { return this->operator()(Return_base_tag(), a, b, c, d); }

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return this->operator()(Return_base_tag(), p, q, r); }

    result_type
    operator()(const Point_3& p, const Direction_3& d) const
    { return this->operator()(Return_base_tag(), p, d); }

    result_type
    operator()(const Point_3& p, const Vector_3& v) const
    { return this->operator()(Return_base_tag(), p, v); }

    result_type
    operator()(const Line_3& l, const Point_3& p) const
    { return this->operator()(Return_base_tag(), l, p); }

    result_type
    operator()(const Ray_3& r, const Point_3& p) const
    { return this->operator()(Return_base_tag(), r, p); }

    result_type
    operator()(const Segment_3& s, const Point_3& p) const
    { return this->operator()(Return_base_tag(), s, p); }

    result_type
    operator()(const Circle_3 & c) const
    { return this->operator()(Return_base_tag(), c); }

    result_type
    operator() ( const typename SK::Polynomial_1_3 &eq )
    { return SphericalFunctors::construct_plane_3<SK>(eq); }

    result_type operator() (const Circular_arc_3 & c) const
    { return c.rep().supporting_plane(); }
  };



  template <class SK>
  class Construct_line_3
  {
    typedef typename SK::Line_arc_3                Line_arc_3;
    typedef typename SK::Line_3                    Line_3;
    typedef typename SK::Point_3 Point_3;
    typedef typename SK::Direction_3 Direction_3;
    typedef typename SK::Vector_3 Vector_3;
    typedef typename SK::Segment_3 Segment_3;
    typedef typename SK::Ray_3 Ray_3;
  public:
    typedef typename SK::Linear_kernel::Construct_line_3 LK_Construct_line_3;
    typedef typename LK_Construct_line_3::result_type result_type;

    result_type
    operator()(Return_base_tag, const Point_3& p, const Point_3& q) const
    { return LK_Construct_line_3()(p, Vector_3(p, q)); }

    result_type
    operator()(Return_base_tag, const Point_3& p, const Direction_3& d) const
    { return operator()(Return_base_tag(), p, Vector_3(d.dx(), d.dy(), d.dz())); }

    result_type
    operator()(Return_base_tag, const Point_3& p, const Vector_3& v) const
    { return LK_Construct_line_3()(p, v); }

    result_type
    operator()(Return_base_tag, const Segment_3& s) const
    { return LK_Construct_line_3()(s.source(), Vector_3(s.source(), s.target())); }

    result_type
    operator()(Return_base_tag, const Ray_3& r) const
    { return LK_Construct_line_3()(r.source(), Vector_3(r.source(), r.second_point())); }


    result_type
    operator()(const Point_3& p, const Point_3& q) const
    { return this->operator()(Return_base_tag(), p, q); }

    result_type
    operator()(const Point_3& p, const Direction_3& d) const
    { return this->operator()(Return_base_tag(), p, d); }

    result_type
    operator()(const Point_3& p, const Vector_3& v) const
    { return this->operator()(Return_base_tag(), p, v); }

    result_type
    operator()(const Segment_3& s) const
    { return this->operator()(Return_base_tag(), s); }

    result_type
    operator()(const Ray_3& r) const
    { return this->operator()(Return_base_tag(), r); }

    const result_type&
    operator() (const Line_arc_3 & a) const
    { return (a.rep().supporting_line()); }

    result_type
    operator() ( const typename SK::Polynomials_for_line_3 &eq )
    { return SphericalFunctors::construct_line_3<SK>(eq); }

  };

  template < class SK >
  class Construct_circle_3
  {
    typedef typename SK::Linear_kernel::Construct_circle_3 Extended;
    typedef typename Extended::result_type forwarded_result_type;
    typedef typename SK::FT FT;
    typedef typename SK::Point_3 Point_3;
    typedef typename SK::Plane_3 Plane_3;
    typedef typename SK::Sphere_3 Sphere_3;
    typedef typename SK::Circle_3 Circle_3;
    typedef typename SK::Vector_3    Vector_3;
    typedef typename SK::Direction_3 Direction_3;
    typedef typename SK::Circular_arc_3            Circular_arc_3;
    typedef typename SK::Kernel_base::Circle_3  RCircle_3;
    typedef typename Circle_3::Rep              Rep;
  public:
    template<typename>
    struct result {
      typedef forwarded_result_type type;
    };

    template<typename F>
    struct result<F(Circular_arc_3)> {
      typedef const forwarded_result_type& type;
    };

    forwarded_result_type
    operator()(const Point_3& p, const FT& sr,
               const Plane_3& plane) const
    { return Extended()(p, sr, plane); }

    forwarded_result_type
    operator() (const Point_3& p, const FT& sr,
                const Vector_3& v) const
    { return Extended()(p, sr, v); }

    forwarded_result_type
    operator() (const Point_3& p, const FT& sr,
                const Direction_3& d) const
    { return Extended()(p, sr, d); }

    forwarded_result_type
    operator() (const Sphere_3& s1, const Sphere_3& s2) const
    { return Extended()(s1, s2); }

    forwarded_result_type
    operator() (const Plane_3& p, const Sphere_3& s) const
    { return Extended()(p, s); }

    forwarded_result_type
    operator() (const Sphere_3& s, const Plane_3& p) const
    { return Extended()(p, s); }

    forwarded_result_type
    operator() (const Plane_3& p, const Sphere_3& s, int a) const
    { return Extended()(p, s, a); }

    forwarded_result_type
    operator() (const Sphere_3& s, const Plane_3& p, int a) const
    { return Extended()(p, s, a); }

    forwarded_result_type
    operator()(        const Point_3& p1, const Point_3& p2, const Point_3& p3) const
    { return Extended()(p1, p2, p3); }

    forwarded_result_type
    operator() ( const typename SK::Polynomials_for_circle_3 &eq )
    { return Rep(construct_circle_3<SK>(eq)); }

    const forwarded_result_type&
    operator() (const Circular_arc_3 & a) const
    { return (a.rep().supporting_circle()); }

  };

  template < class SK >
  class Construct_line_arc_3
  {

    typedef typename SK::Line_3                    Line_3;
    typedef typename SK::Point_3                   Point_3;
    typedef typename SK::Segment_3                 Segment_3;
    typedef typename SK::Sphere_3                  Sphere_3;
    typedef typename SK::Plane_3                   Plane_3;
    typedef typename SK::Circular_arc_point_3      Circular_arc_point_3;
    typedef typename SK::Line_arc_3                Line_arc_3;
    typedef typename SK::Kernel_base::Line_arc_3   RLine_arc_3;
    typedef typename Line_arc_3::Rep               Rep;

  public:
    typedef Line_arc_3   result_type;

    result_type
    operator()(void) const
    { return Rep(); }

    result_type
    operator()(const Line_3 &l,
               const Circular_arc_point_3 &s,
               const Circular_arc_point_3 &t) const
    { return Rep(l,s,t); }

    result_type
    operator()(const Point_3 &s,
                     const Point_3 &t) const
    { return Rep(s,t); }

    result_type
    operator()(const Segment_3 &s) const
    { return Rep(s); }

    // Not Documented
    result_type
    operator()(const Line_3 &l,
               const Sphere_3 &s,
               bool less_xyz_first = true) const
    { return Rep(l,s,less_xyz_first); }

    // Not Documented
    result_type
    operator()(const Sphere_3 &s,
               const Line_3 &l,
               bool less_xyz_first = true) const
    { return Rep(l,s,less_xyz_first); }

    // Not Documented
    result_type
    operator()(const Line_3 &l,
               const Sphere_3 &s1, bool less_xyz_s1,
               const Sphere_3 &s2, bool less_xyz_s2) const
    { return Rep(l,s1,less_xyz_s1,
                   s2,less_xyz_s2); }

    // Not Documented
    result_type
    operator()(const Sphere_3 &s1, bool less_xyz_s1,
               const Sphere_3 &s2, bool less_xyz_s2,
               const Line_3 &l) const
    { return Rep(l,s1,less_xyz_s1,
                   s2,less_xyz_s2); }

    // Not Documented
    result_type
    operator()(const Line_3 &l,
               const Plane_3 &p1,
               const Plane_3 &p2) const
    { return Rep(l,p1,p2); }

    // Not Documented
    result_type
    operator()(const Plane_3 &p1,
               const Plane_3 &p2,
               const Line_3 &l) const
    { return Rep(l,p1,p2); }

  };

  template < class SK >
  class Construct_circular_arc_3
  {

    typedef typename SK::Line_3                        Line_3;
    typedef typename SK::Point_3                       Point_3;
    typedef typename SK::Segment_3                     Segment_3;
    typedef typename SK::Sphere_3                      Sphere_3;
    typedef typename SK::Plane_3                       Plane_3;
    typedef typename SK::Line_arc_3                    Line_arc_3;

    typedef typename SK::Circular_arc_point_3          Circular_arc_point_3;
    typedef typename SK::Circle_3                      Circle_3;
    typedef typename SK::Circular_arc_3                Circular_arc_3;
    typedef typename SK::Kernel_base::Circular_arc_3   RCircular_arc_3;
    typedef typename Circular_arc_3::Rep               Rep;
  public:
    typedef Circular_arc_3   result_type;

    result_type
    operator()(void) const
    { return Rep(); }

    result_type
    operator()(const Circle_3 &c) const
    { return Rep(c); }

    result_type
    operator()(const Circle_3 &c,const Circular_arc_point_3& pt) const
    { return Rep(c,pt); }

    result_type
    operator()(const Circle_3 &l,
               const Circular_arc_point_3 &s,
               const Circular_arc_point_3 &t) const
    { return Rep(l,s,t); }

    // Not Documented
    result_type
    operator()(const Circle_3 &c,
               const Sphere_3 &s1, bool less_xyz_s1,
               const Sphere_3 &s2, bool less_xyz_s2) const
    { return Rep(c,s1,less_xyz_s1,s2,less_xyz_s2); }

    // Not Documented
    result_type
    operator()(const Sphere_3 &s1, bool less_xyz_s1,
               const Sphere_3 &s2, bool less_xyz_s2,
               const Circle_3 &c) const
    { return Rep(c,s1,less_xyz_s1,s2,less_xyz_s2); }

    // Not Documented
    result_type
    operator()(const Circle_3 &c,
               const Plane_3 &p1, bool less_xyz_p1,
               const Plane_3 &p2, bool less_xyz_p2) const
    { return Rep(c,p1,less_xyz_p1,p2,less_xyz_p2); }

    // Not Documented
    result_type
    operator()(const Plane_3 &p1, bool less_xyz_p1,
               const Plane_3 &p2, bool less_xyz_p2,
               const Circle_3 &c) const
    { return Rep(c,p1,less_xyz_p1,p2,less_xyz_p2); }

    result_type
    operator()(const Point_3 &begin,
               const Point_3 &middle,
               const Point_3 &end) const
    { return Rep(begin,middle,end); }

  };

  template <class SK>
  class Construct_circular_min_vertex_3
  {
    typedef typename SK::Line_arc_3                Line_arc_3;
    typedef typename SK::Circular_arc_point_3      Circular_arc_point_3;

  public:

    typedef const Circular_arc_point_3& result_type;

    result_type operator() (const Line_arc_3 & a) const
    { return (a.rep().lower_xyz_extremity()); }

  };

  template <class SK>
  class Construct_circular_max_vertex_3
  {
    typedef typename SK::Line_arc_3                Line_arc_3;
    typedef typename SK::Circular_arc_point_3      Circular_arc_point_3;

  public:

    typedef const Circular_arc_point_3& result_type;

    result_type operator() (const Line_arc_3 & a) const
    { return (a.rep().higher_xyz_extremity()); }

  };

  template <class SK>
  class Construct_circular_source_vertex_3
  {
    typedef typename SK::Line_arc_3                Line_arc_3;
    typedef typename SK::Circular_arc_3            Circular_arc_3;
    typedef typename SK::Circular_arc_point_3      Circular_arc_point_3;

  public:

    typedef const Circular_arc_point_3& result_type;

    result_type operator() (const Line_arc_3 & a) const
    { return (a.rep().source()); }

    result_type operator() (const Circular_arc_3 & a) const
    { return (a.rep().source()); }

  };

  template <class SK>
  class Construct_circular_target_vertex_3
  {
    typedef typename SK::Line_arc_3                Line_arc_3;
    typedef typename SK::Circular_arc_point_3      Circular_arc_point_3;
    typedef typename SK::Circular_arc_3            Circular_arc_3;

  public:

    typedef const Circular_arc_point_3& result_type;

    result_type operator() (const Line_arc_3 & a) const
    { return (a.rep().target()); }

    result_type operator() (const Circular_arc_3 & a) const
    { return (a.rep().target()); }

  };

  template < class SK >
  class Has_on_3
    : public SK::Linear_kernel::Has_on_3
  {
    typedef typename SK::Point_3                 Point_3;
    typedef typename SK::Sphere_3                Sphere_3;
    typedef typename SK::Plane_3                 Plane_3;
    typedef typename SK::Line_3                  Line_3;
    typedef typename SK::Segment_3                   Segment_3;
    typedef typename SK::Ray_3                   Ray_3;
    typedef typename SK::Triangle_3              Triangle_3;
    typedef typename SK::Line_arc_3              Line_arc_3;
    typedef typename SK::Circular_arc_point_3    Circular_arc_point_3;
    typedef typename SK::Circular_arc_3          Circular_arc_3;
    typedef typename SK::Circle_3                Circle_3;


  public:
    typedef typename SK::Linear_kernel::Has_on_3::result_type result_type;

    using SK::Linear_kernel::Has_on_3::operator();

    result_type
    operator()(const Sphere_3 &a, const Circular_arc_point_3 &p) const
    { return SphericalFunctors::has_on<SK>(a, p); }

    result_type
    operator()(const Plane_3 &a, const Circular_arc_point_3 &p) const
    { return SphericalFunctors::has_on<SK>(a, p); }

    result_type
    operator()(const Line_3 &a, const Circular_arc_point_3 &p) const
    { return SphericalFunctors::has_on<SK>(a, p); }

    result_type
    operator()(const Circle_3 &a, const Circular_arc_point_3 &p) const
    { return SphericalFunctors::has_on<SK>(a, p); }

    result_type
    operator()(const Line_arc_3 &a, const Circular_arc_point_3 &p,
               const bool already_know_point_on_line = false) const
    { return SphericalFunctors::has_on<SK>(a, p, already_know_point_on_line); }

    result_type
    operator()(const Line_arc_3 &a, const Point_3 &p,
               const bool already_know_point_on_line = false) const
    { return SphericalFunctors::has_on<SK>(a, p, already_know_point_on_line); }

    result_type
    operator()(const Plane_3 &p, const Line_arc_3 &a) const
    { return SphericalFunctors::has_on<SK>(p, a); }

    result_type
    operator()(const Line_3 &a, const Line_arc_3 &p) const
    { return SphericalFunctors::has_on<SK>(a, p); }

    result_type
    operator()(const Circular_arc_3 &a, const Point_3 &p,
               const bool has_on_supporting_circle = false) const
    { return SphericalFunctors::has_on<SK>(a, p, has_on_supporting_circle); }

    result_type
    operator()(const Circular_arc_3 &a, const Circular_arc_point_3 &p,
               const bool has_on_supporting_circle = false) const
    { return SphericalFunctors::has_on<SK>(a, p, has_on_supporting_circle); }

    result_type
    operator()(const Sphere_3 &a, const Circular_arc_3 &p) const
    { return SphericalFunctors::has_on<SK>(a, p); }

    result_type
    operator()(const Plane_3 &a, const Circular_arc_3 &p) const
    { return SphericalFunctors::has_on<SK>(a, p); }

    result_type
    operator()(const Circle_3 &a, const Circular_arc_3 &p) const
    { return SphericalFunctors::has_on<SK>(a, p); }

    result_type
    operator()(const Circular_arc_3 &p, const Circle_3 &a) const
    { return SphericalFunctors::has_on<SK>(p, a); }

  };

  template < class SK >
  class Do_intersect_3
    : public SK::Linear_kernel::Do_intersect_3
  {

    typedef typename SK::Sphere_3                 Sphere_3;
    typedef typename SK::Line_3                   Line_3;
    typedef typename SK::Line_arc_3               Line_arc_3;
    typedef typename SK::Circular_arc_3           Circular_arc_3;
    typedef typename SK::Plane_3                  Plane_3;
    typedef typename SK::Circle_3                 Circle_3;
    typedef typename SK::Circle_3                 Circular_arc_point_3;

  public:
    typedef typename SK::Linear_kernel::Do_intersect_3::result_type result_type;

    using SK::Linear_kernel::Do_intersect_3::operator();

#define CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(A,B)            \
    result_type                                                         \
    operator()(const A & c1, const B & c2) const                        \
    { std::vector< typename SK3_Intersection_traits<SK, A, B>::type > res; \
      typename SK::Intersect_3()(c1,c2,std::back_inserter(res));        \
      return !res.empty(); }

#define CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_3(A,B,C)          \
    result_type                                                         \
    operator()(const A & c1, const B & c2, const C & c3) const          \
    { std::vector< typename SK3_Intersection_traits<SK, A, B, C>::type > res; \
      typename SK::Intersect_3()(c1,c2,c3,std::back_inserter(res));     \
      return !res.empty(); }

        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Sphere_3, Line_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_3, Sphere_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_3(Sphere_3, Sphere_3, Sphere_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_3(Sphere_3, Sphere_3, Plane_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_3(Plane_3, Sphere_3, Sphere_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_3(Plane_3, Plane_3, Sphere_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_3(Sphere_3, Plane_3, Plane_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circle_3, Plane_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Plane_3, Circle_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circle_3, Sphere_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Sphere_3, Circle_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circle_3, Circle_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circle_3, Line_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_3, Circle_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_arc_3, Line_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_3, Line_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_arc_3, Line_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circle_3, Line_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_arc_3, Circle_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Sphere_3, Line_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_arc_3, Sphere_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Plane_3, Line_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_arc_3, Plane_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circular_arc_3, Circular_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_3, Circular_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circular_arc_3, Line_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circle_3, Circular_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circular_arc_3, Circle_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Sphere_3, Circular_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circular_arc_3, Sphere_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Plane_3, Circular_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circular_arc_3, Plane_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Circular_arc_3, Line_arc_3)
        CGAL_SPHERICAL_KERNEL_MACRO_DO_INTERSECTION_3_2(Line_arc_3, Circular_arc_3)

  };

  template < class SK >
  class Intersect_3
  //The inheritance is commented as for some reason this does not work when
  //using the Lazy_kernel as linear kernel.
  //  : public SK::Linear_kernel::Intersect_3
  {
    typedef typename SK::Sphere_3                 Sphere_3;
    typedef typename SK::Line_3                   Line_3;
    typedef typename SK::Line_arc_3               Line_arc_3;
    typedef typename SK::Circular_arc_3           Circular_arc_3;
    typedef typename SK::Plane_3                  Plane_3;
    typedef typename SK::Point_3                  Point_3;
    typedef typename SK::Circle_3                 Circle_3;
    typedef typename SK::Circular_arc_point_3     Circular_arc_point_3;

  public:

    template <typename>
    struct result;

    // the binary overload always goes to Linear::Intersect_3
    template <typename F, typename A, typename B>
    struct result<F(A, B)>
    { typedef typename Intersection_traits<SK, A, B>::result_type type; };

    // This one is only for the spherical kernel, O is an output iterator
    template <typename F, typename A, typename B, typename OutputIterator>
    struct result<F(A, B, OutputIterator)>
    { typedef OutputIterator type;};

    // there is no quaternary form in the linear Kernel
    template <typename F, typename A, typename B, typename C, typename OutputIterator>
    struct result<F(A, B, C, OutputIterator)>
    { typedef OutputIterator type; };

    //only ternary from the linear kernel
    template<typename F>
    struct result<F(Plane_3, Plane_3, Plane_3)> {
      typedef boost::optional<
        boost::variant< Point_3,
                        Line_3,
                        Plane_3 > > type;
    };

    //using SK::Linear_kernel::Intersect_3::operator();

    typedef typename SK::Linear_kernel::Intersect_3 Intersect_linear_3;

    template<class A, class B>
    typename Intersection_traits<SK, A, B>::result_type
    operator()(const A& a, const B& b) const{
      return Intersect_linear_3()(a,b);
    }

    typename result<Intersect_linear_3(Plane_3, Plane_3, Plane_3)>::type
    operator()(const Plane_3& p, const Plane_3& q, const Plane_3& r) const
    {
      return Intersect_linear_3()(p, q, r);
    }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s, const Line_3 & l,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (s,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_3 & l,const Sphere_3 & s,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (s,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s1, const Sphere_3 & s2,
               const Sphere_3 & s3, OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (s1,s2,s3,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s1, const Sphere_3 & s2,
               const Plane_3 & p, OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (p,s1,s2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Plane_3 & p, const Sphere_3 & s1,
               const Sphere_3 & s2, OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (p,s1,s2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Plane_3 & p1, const Plane_3 & p2,
               const Sphere_3 & s, OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (p1,p2,s,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s, const Plane_3 & p1,
               const Plane_3 & p2, OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (p1,p2,s,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c, const Plane_3 & p,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,p,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Plane_3 & p, const Circle_3 & c,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,p,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c, const Sphere_3 & s,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,s,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s, const Circle_3 & c,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,s,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c1, const Circle_3 & c2,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c1,c2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c, const Line_3 & l,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_3 & l, const Circle_3 & c,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_3 & l1, const Line_arc_3 & l2,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (l1,l2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_3 & l, const Line_arc_3 & la,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (l,la,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_3 & la, const Line_3 & l,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (l,la,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c, const Line_arc_3 & l,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_3 & l, const Circle_3 & c,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s, const Line_arc_3 & l,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (s,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_3 & l,const Sphere_3 & s,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (s,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Plane_3 & s, const Line_arc_3 & l,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (s,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_3 & l,const Plane_3 & s,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (s,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_3 & c1, const Circular_arc_3 & c2,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c1,c2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_3 & l, const Circular_arc_3 & ca,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (l,ca,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_3 & ca, const Line_3 & l,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (l,ca,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c, const Circular_arc_3 & ca,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,ca,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_3 & ca, const Circle_3 & c,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (c,ca,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s, const Circular_arc_3 & ca,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (s,ca,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_3 & ca,const Sphere_3 & s,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (s,ca,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Plane_3 & p, const Circular_arc_3 & ca,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (p,ca,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_3 & ca, const Plane_3 & p,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (p,ca,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_3 & ca, const Line_arc_3 & la,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (ca,la,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_3 & la, const Circular_arc_3 & ca,
               OutputIterator res) const
    { return SphericalFunctors::intersect_3<SK> (ca,la,res); }

  };

// If 2 line_arc have the same supporting line
// if they intersect only at a point
// even in this case we consider that the 2 line_arc overlap
  template < class SK >
  class Do_overlap_3
  {
    typedef typename SK::Line_arc_3     Line_arc_3;
    typedef typename SK::Circular_arc_3 Circular_arc_3;

  public:
    typedef bool         result_type;

    result_type
    operator() (const Line_arc_3 &l1, const Line_arc_3 &l2,
                const bool known_equal_supporting_line = false) const
    { return SphericalFunctors::do_overlap<SK>(l1, l2, known_equal_supporting_line); }

    result_type
    operator() (const Circular_arc_3 &c1, const Circular_arc_3 &c2,
                const bool known_equal_supporting_circle = false) const
    { return SphericalFunctors::do_overlap<SK>(c1, c2, known_equal_supporting_circle); }

  };

  template < class SK >
  class Split_3
  {
    typedef typename SK::Circular_arc_point_3    Circular_arc_point_3;
    typedef typename SK::Line_arc_3              Line_arc_3;
    typedef typename SK::Circular_arc_3          Circular_arc_3;

  public:

    typedef void         result_type;

    result_type
    operator()(const Line_arc_3 &l,
               const Circular_arc_point_3 &p,
               Line_arc_3 &ca1, Line_arc_3 &ca2) const
    { return SphericalFunctors::split<SK>(l, p, ca1, ca2); }

    result_type
    operator()(const Circular_arc_3 &c,
               const Circular_arc_point_3 &p,
               Circular_arc_3 &ca1, Circular_arc_3 &ca2) const
    { return SphericalFunctors::split<SK>(c, p, ca1, ca2); }

  };

  template <class SK>
  class Construct_bbox_3
    : public SK::Linear_kernel::Construct_bbox_3
  {
    typedef typename SK::Circular_arc_point_3      Circular_arc_point_3;
    typedef typename SK::Circular_arc_3            Circular_arc_3;
    typedef typename SK::Point_3                   Point_3;
    typedef typename SK::Segment_3                 Segment_3;
    typedef typename SK::Circle_3                  Circle_3;
    typedef typename SK::Triangle_3                Triangle_3;
    typedef typename SK::Tetrahedron_3             Tetrahedron_3;
    typedef typename SK::Sphere_3                  Sphere_3;
    typedef typename SK::Iso_cuboid_3              Iso_cuboid_3;
    typedef typename SK::Line_arc_3                Line_arc_3;

  public:

    typedef typename SK::Linear_kernel::Construct_bbox_3::result_type result_type;

    using SK::Linear_kernel::Construct_bbox_3::operator();

    result_type operator() (const Circular_arc_point_3 & c) const
    { return c.rep().bbox(); }

    result_type operator() (const Line_arc_3 & l) const
    { return l.rep().bbox(); }

    result_type operator() (const Circular_arc_3 & c) const
    { return c.rep().bbox(); }

  };

  template <class SK>
  class Compute_approximate_squared_length_3
    : public SK::Linear_kernel::Compute_approximate_squared_length_3
  {
    typedef typename SK::Circle_3                  Circle_3;
    typedef typename SK::Circular_arc_3            Circular_arc_3;
    typedef typename SK::FT                        FT;

  public:

    typedef double result_type;
    using SK::Linear_kernel::Compute_approximate_squared_length_3::operator();

    result_type operator() (const Circular_arc_3 & c) const
    { return c.rep().approximate_squared_length(); }

  };

  template <class SK>
  class Compute_approximate_angle_3
    : public SK::Linear_kernel::Compute_approximate_angle_3
  {
    typedef typename SK::Point_3                   Point_3;
    typedef typename SK::Vector_3                  Vector_3;
    typedef typename SK::Circular_arc_3            Circular_arc_3;
    typedef typename SK::FT                        FT;

  public:

    typedef double result_type;

    using SK::Linear_kernel::Compute_approximate_angle_3::operator();

    result_type operator() (const Circular_arc_3 & c) const
    { return c.rep().approximate_angle(); }

  };

  template <class SK>
  class Bounded_side_3
    : public SK::Linear_kernel::Bounded_side_3
  {
    typedef typename SK::Sphere_3              Sphere_3;
    typedef typename SK::Circle_3              Circle_3;
    typedef typename SK::Circular_arc_point_3  Circular_arc_point_3;
    typedef typename SK::Point_3               Point_3;

  public:
    typedef typename SK::Linear_kernel::Bounded_side_3::result_type    result_type;

    using SK::Linear_kernel::Bounded_side_3::operator();

    result_type
    operator()( const Sphere_3& s, const Circular_arc_point_3& p) const
    { return SphericalFunctors::bounded_side<SK>(s,p); }

    result_type
    operator()( const Circle_3& c, const Circular_arc_point_3& p) const
    { return SphericalFunctors::bounded_side<SK>(c,p); }

    // We can maybe optimize it doing the operator() for point_3 too

  };

  template <class SK>
  class Has_on_bounded_side_3
    : public SK::Linear_kernel::Has_on_bounded_side_3
  {
    typedef typename SK::Sphere_3              Sphere_3;
    typedef typename SK::Circle_3              Circle_3;
    typedef typename SK::Circular_arc_point_3  Circular_arc_point_3;
    typedef typename SK::Point_3               Point_3;

  public:
    typedef typename SK::Linear_kernel::Has_on_bounded_side_3::result_type    result_type;

    using SK::Linear_kernel::Has_on_bounded_side_3::operator();

    result_type
    operator()( const Sphere_3& s, const Circular_arc_point_3& p) const
    { return SK().bounded_side_3_object()(s,p) == ON_BOUNDED_SIDE; }

    result_type
    operator()( const Circle_3& c, const Circular_arc_point_3& p) const
    { return SK().bounded_side_3_object()(c,p) == ON_BOUNDED_SIDE; }

    // We can maybe optimize it doing the operator() for point_3 too

  };

  template <class SK>
  class Has_on_unbounded_side_3
    : public SK::Linear_kernel::Has_on_unbounded_side_3
  {
    typedef typename SK::Sphere_3              Sphere_3;
    typedef typename SK::Circle_3              Circle_3;
    typedef typename SK::Circular_arc_point_3  Circular_arc_point_3;
    typedef typename SK::Point_3               Point_3;

  public:
    typedef typename SK::Linear_kernel::Has_on_unbounded_side_3::result_type    result_type;

    using SK::Linear_kernel::Has_on_unbounded_side_3::operator();

    result_type
    operator()( const Sphere_3& s, const Circular_arc_point_3& p) const
    { return SK().bounded_side_3_object()(s,p) == ON_UNBOUNDED_SIDE; }

    result_type
    operator()( const Circle_3& c, const Circular_arc_point_3& p) const
    { return SK().bounded_side_3_object()(c,p) == ON_UNBOUNDED_SIDE; }

    // We can maybe optimize it doing the operator() for point_3 too

  };

  template <class SK>
  class Is_theta_monotone_3{
    typename SK::Sphere_3 sphere_;
  public:
    typedef bool result_type;

    Is_theta_monotone_3(const typename SK::Sphere_3& sphere):sphere_(sphere){}

    result_type
    operator()(const typename SK::Circular_arc_3& arc) const {
      return SphericalFunctors::is_theta_monotone_3<SK>(arc,sphere_);
    }
  };

  template < class SK >
  class Compare_theta_3{
    typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
    typedef typename SK::Vector_3 Vector_3;

    typename SK::Sphere_3 sphere_;

  public:
    typedef  CGAL::Comparison_result result_type;

    Compare_theta_3(const typename SK::Sphere_3& sphere):sphere_(sphere){}

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return SphericalFunctors::compare_theta_of_pts<SK>(p0, p1,sphere_); }

    result_type
    operator() (const Circular_arc_point_3 &p,
                const Vector_3 &v) const
    { return SphericalFunctors::compare_theta_pt_vector<SK>(p,v,sphere_); }

    result_type
    operator() (const Vector_3 &m1,
                const Vector_3 &m2) const
    { return SphericalFunctors::compare_theta_vectors<SK>(m1,m2); }

    result_type
    operator() (const Vector_3 &v,const Circular_arc_point_3 &p0) const
    { return CGAL::opposite( SphericalFunctors::compare_theta_pt_vector<SK>(p0, v,sphere_) ); }
  };

  template < class SK >
  class Compare_theta_z_3{
    typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
    typename SK::Sphere_3 sphere_;

  public:
    typedef  CGAL::Comparison_result result_type;

    Compare_theta_z_3(const typename SK::Sphere_3& sphere):sphere_(sphere){}

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Circular_arc_point_3 &p1,bool decreasing_z=false) const
    { return SphericalFunctors::compare_theta_z<SK>(p0, p1,sphere_,decreasing_z); }

  };

  template < class SK >
  class Make_theta_monotone_3{
    typename SK::Sphere_3 sphere_;

  public:
    Make_theta_monotone_3(const typename SK::Sphere_3& sphere):sphere_(sphere){}

    template <class OutputIterator>
    OutputIterator
    operator() (const typename SK::Circle_3 &circle,OutputIterator out_it) const
    { return SphericalFunctors::make_circle_theta_monotone<SK>(circle,sphere_,out_it); }

    template <class OutputIterator>
    OutputIterator
    operator() (const typename SK::Circular_arc_3 &arc,OutputIterator out_it) const
    { return SphericalFunctors::make_circular_arc_theta_monotone<SK>(arc,sphere_,out_it); }

  };

  template <class SK>
  class Compare_z_to_right_3{
    typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
    typedef typename SK::Circular_arc_3 Circular_arc_3;
    typename SK::Sphere_3 sphere_;

  public:
    typedef  CGAL::Comparison_result result_type;
    Compare_z_to_right_3(const typename SK::Sphere_3& sphere):sphere_(sphere){}

    result_type
    operator()(const Circular_arc_3& arc1,const Circular_arc_3& arc2,const Circular_arc_point_3& pt,bool do_to_the_left=false){
      CGAL_kernel_precondition(SK().has_on_3_object()(sphere_,arc1));
      CGAL_kernel_precondition(SK().has_on_3_object()(sphere_,arc2));
      CGAL_kernel_precondition(SphericalFunctors::is_theta_monotone_3<SK>(arc1,sphere_));
      CGAL_kernel_precondition(SphericalFunctors::is_theta_monotone_3<SK>(arc2,sphere_));
      CGAL_kernel_precondition(SK().has_on_3_object()(arc1,pt));
      CGAL_kernel_precondition(SK().has_on_3_object()(arc2,pt));
      CGAL_kernel_precondition(classify_circle_3<SK>(arc1.supporting_circle(),sphere_)!=NORMAL || pt.z()!=extremal_points_z_coordinate<SK>(arc1.supporting_circle(),sphere_));
      CGAL_kernel_precondition(classify_circle_3<SK>(arc2.supporting_circle(),sphere_)!=NORMAL || pt.z()!=extremal_points_z_coordinate<SK>(arc2.supporting_circle(),sphere_));

      if (pt.y() == sphere_.center().y() && pt.x() > sphere_.center().x()) //case theta = 0
        return
          Compare_to_right_of_arcs<SK,Trait_for_cmp_tgt_theta_0<SK> >( Trait_for_cmp_tgt_theta_0<SK>(pt.coordinates(),sphere_),sphere_ )
            (arc1,arc2,do_to_the_left);
      else //general case
        return
          Compare_to_right_of_arcs<SK,Trait_for_cmp_tgt<SK> >( Trait_for_cmp_tgt<SK>(pt.coordinates(),sphere_),sphere_ )
            (arc1,arc2,do_to_the_left);
    }
  };

  template <class SK>
  class Compare_z_at_theta_3{
     typename SK::Sphere_3 sphere_;

  public:
    typedef  CGAL::Comparison_result result_type;
    Compare_z_at_theta_3(const typename SK::Sphere_3& sphere):sphere_(sphere){}

    result_type
    operator()( const typename SK::Circular_arc_3& arc1,
                const typename SK::Circular_arc_3& arc2,
                const typename SK::Vector_3& m) const
    {
      return SphericalFunctors::compare_z_at_theta_arcs<SK>(arc1,arc2,m,sphere_);
    }

    result_type
    operator()(const typename SK::Circular_arc_point_3& point,
               const typename SK::Circular_arc_3& arc) const
    {
      return SphericalFunctors::compare_z_at_theta_pt_arc<SK>(point,arc,sphere_);
    }
  };


} // namespace SphericalFunctors

} // namespace CGAL

#endif // CGAL_SPHERICAL_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_SPHERE_H

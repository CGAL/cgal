// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
//           Damien Leroy
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCLE_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCLE_3_H

namespace CGAL {
  namespace SphericalFunctors {
            
    template < class SK >
    typename SK::Circle_3
    construct_circle_3(const typename SK::Polynomials_for_circle_3 &eq)
    {
      typedef typename SK::Point_3 Point_3;
      typedef typename SK::Plane_3 Plane_3;
      typedef typename SK::Circle_3 Circle_3;
      typedef typename SK::Sphere_3 Sphere_3;
      typedef typename SK::FT FT;
      Sphere_3 s = construct_sphere_3<SK>(eq.first);
      Plane_3 p = construct_plane_3<SK>(eq.second);
      const FT d2 = CGAL::square(p.a()*s.center().x() + 
                                 p.b()*s.center().y() + 
                                 p.c()*s.center().z() + p.d()) /
       (CGAL::square(p.a()) + CGAL::square(p.b()) + CGAL::square(p.c()));
      // We do not accept circles with radius 0 (should we?)
      CGAL_kernel_precondition(d2 <  s.squared_radius());
      // d2 < s.squared_radius()
      Point_3 center = p.projection(s.center());
      return Circle_3(center,s.squared_radius() - d2,p);
    }

    template< class SK>
    bool
    equal( const typename SK::Circle_3 &p1,
           const typename SK::Circle_3 &p2)
    {
      return p1.rep() == p2.rep();
    }

    template <class SK>
    typename SK::Circular_arc_point_3
    x_extremal_point(const typename SK::Circle_3 & c, bool i)
    {
      typedef typename SK::Algebraic_kernel   AK;
      return AK().x_critical_points_object()(typename SK::Get_equation()(c),i);
    }

    template <class SK,class OutputIterator>
    OutputIterator
    x_extremal_points(const typename SK::Circle_3 & c, OutputIterator res)
    {
      typedef typename SK::Algebraic_kernel   AK;
      return AK().x_critical_points_object()(typename SK::Get_equation()(c),res);
    }

    template <class SK>
    typename SK::Circular_arc_point_3
    y_extremal_point(const typename SK::Circle_3 & c, bool i)
    {
      typedef typename SK::Algebraic_kernel   AK;
      return AK().y_critical_points_object()(typename SK::Get_equation()(c),i);
    }

    template <class SK,class OutputIterator>
    OutputIterator
    y_extremal_points(const typename SK::Circle_3 & c, OutputIterator res)
    {
      typedef typename SK::Algebraic_kernel   AK;
      return AK().y_critical_points_object()(typename SK::Get_equation()(c),res);
    }

    template <class SK>
    typename SK::Circular_arc_point_3
    z_extremal_point(const typename SK::Circle_3 & c, bool i)
    {
      typedef typename SK::Algebraic_kernel   AK;
      return AK().z_critical_points_object()(typename SK::Get_equation()(c),i);
    }

    template <class SK,class OutputIterator>
    OutputIterator
    z_extremal_points(const typename SK::Circle_3 & c, OutputIterator res)
    {
      typedef typename SK::Algebraic_kernel   AK;
      return AK().z_critical_points_object()(typename SK::Get_equation()(c),res);
    }

  }//SphericalFunctors
}//CGAL



#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCLE_3_H

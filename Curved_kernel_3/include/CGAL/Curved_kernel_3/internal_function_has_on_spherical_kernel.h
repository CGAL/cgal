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

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_HAS_ON_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_HAS_ON_3_H

namespace CGAL {
  namespace SphericalFunctors {

    template <class SK>
    inline
    bool
    has_on(const typename SK::Sphere_3 &a, 
           const typename SK::Point_3 &p)
    { 
      return a.has_on_boundary(p);
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Sphere_3 &a, 
           const typename SK::Circular_arc_point_3 &p)
    { 
      typedef typename SK::AK AK;
      typedef typename SK::Polynomial_for_spheres_2_3 Equation;
      Equation equation = get_equation<SK>(a);
      return (AK().sign_at_object()(equation,p.coordinates()) == ZERO);
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Plane_3 &a, 
           const typename SK::Point_3 &p)
    { 
      return a.has_on(p);
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Plane_3 &a, 
           const typename SK::Circular_arc_point_3 &p)
    { 
      typedef typename SK::AK AK;
      typedef typename SK::Polynomial_1_3 Equation;
      Equation equation = get_equation<SK>(a);
      return (AK().sign_at_object()(equation,p.coordinates()) == ZERO);
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Line_3 &a,
           const typename SK::Point_3 &p)
    { 
      return a.has_on(p);
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Line_3 &a, 
           const typename SK::Circular_arc_point_3 &p)
    { 
      typedef typename SK::AK AK;
      typedef typename SK::Polynomials_for_line_3 Equation;
      Equation equation = get_equation<SK>(a);
      return p.coordinates().is_on_line(equation);
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Circle_3 &a, 
           const typename SK::Point_3 &p)
    { 
      return has_on<SK>(a.diametral_sphere(),p) &&
             has_on<SK>(a.supporting_plane(),p);
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Circle_3 &a, 
           const typename SK::Circular_arc_point_3 &p)
    { 
      return has_on<SK>(a.diametral_sphere(),p) &&
             has_on<SK>(a.supporting_plane(),p);
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Sphere_3 &a, 
           const typename SK::Circle_3 &p)
    { 
      typedef typename SK::Point_3 Point_3;
      typedef typename SK::FT FT;
      Point_3 proj = p.supporting_plane().projection(a.center());
      if(!(proj == p.center())) return false;
      const FT d2 = CGAL::square(a.center().x() - p.center().x()) +
                    CGAL::square(a.center().y() - p.center().y()) +
                    CGAL::square(a.center().z() - p.center().z());
      return ((a.squared_radius() - d2) == p.squared_radius());
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Plane_3 &a, 
           const typename SK::Line_3 &p)
    { 
      return a.has_on(p);
    }

    template <class SK>
    inline
    bool
    has_on(const typename SK::Plane_3 &a, 
           const typename SK::Circle_3 &p)
    { 
      return non_oriented_equal<SK>(a,p.supporting_plane());
    }

  }//SphericalFunctors
}//CGAL



#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_HAS_ON_3_H

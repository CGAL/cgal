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



#ifndef CGAL_SPHERICAL_KERNEL_GET_EQUATION_OBJECT_3_H
#define CGAL_SPHERICAL_KERNEL_GET_EQUATION_OBJECT_3_H

#include <CGAL/Curved_kernel_3/internal_functions_on_sphere_3.h>
#include <CGAL/Curved_kernel_3/internal_functions_on_line_3.h>
#include <CGAL/Curved_kernel_3/internal_functions_on_plane_3.h>
#include <CGAL/Curved_kernel_3/internal_functions_on_circle_3.h>
// to be removed when CGAL::Kernel has a Get_equation

namespace CGAL {
  namespace SphericalFunctors {
    
  template < class SK >
  class Get_equation //: public LinearFunctors::Get_equation<SK>
  {
    public:

    typedef void result_type; // should we keep this?

    typedef typename SK::Polynomial_for_spheres_2_3 result_type_for_sphere;
    typedef typename SK::Polynomial_1_3 result_type_for_plane;
    typedef typename SK::Polynomials_for_line_3 result_type_for_line;
    typedef typename SK::Polynomials_for_circle_3 result_type_for_circle;
    //using LinearFunctors::Get_equation<SK>::operator();


    result_type_for_sphere
    operator() ( const typename SK::Sphere_3 & s )
    {
      return get_equation<SK>(s);
    }

    result_type_for_plane
    operator() ( const typename SK::Plane_3 & p )
    {
      return get_equation<SK>(p);
    }

    result_type_for_line
    operator() ( const typename SK::Line_3 & l )
    {
      return get_equation<SK>(l);
    }

    result_type_for_circle
    operator() ( const typename SK::Circle_3 & c )
    {
      return get_equation<SK>(c);
    }

  };
    
  } // namespace SphericalFunctors
} // namespace CGAL

#endif // CGAL_SPHERICAL_KERNEL_GET_EQUATION_OBJECT_3_H

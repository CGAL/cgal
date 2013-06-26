// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_GET_EQUATION_OBJECT_3_H
#define CGAL_SPHERICAL_KERNEL_GET_EQUATION_OBJECT_3_H

#include <CGAL/Circular_kernel_3/internal_functions_on_sphere_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_line_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_plane_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_circle_3.h>
// to be removed when CGAL::Kernel has a Get_equation

namespace CGAL {
  namespace SphericalFunctors {
    
  template < class SK >
  class Get_equation //: public LinearFunctors::Get_equation<SK>
  {
    public:


    typedef typename SK::Polynomial_for_spheres_2_3 result_type_for_sphere;
    typedef typename SK::Polynomial_1_3 result_type_for_plane;
    typedef typename SK::Polynomials_for_line_3 result_type_for_line;
    typedef typename SK::Polynomials_for_circle_3 result_type_for_circle;
    //using LinearFunctors::Get_equation<SK>::operator();

    template <typename>
    struct result;

    template <typename F>
    struct result<F(typename SK::Sphere_3)>
    {
      typedef result_type_for_sphere type;
    };

    template <typename F>
    struct result<F(typename SK::Plane_3)>
    {
      typedef result_type_for_plane type;
    };

    template <typename F>
    struct result<F(typename SK::Line_3)>
    {
      typedef result_type_for_line type;
    };

    template <typename F>
    struct result<F(typename SK::Circle_3)>
    {
      typedef result_type_for_circle type;
    };

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

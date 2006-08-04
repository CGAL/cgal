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

#ifndef CGAL_SPHERICAL_KERNEL_3_H
#define CGAL_SPHERICAL_KERNEL_3_H

#include <CGAL/Curved_kernel_3/Circular_arc_point_3.h>
#include <CGAL/Circular_arc_point_3.h>

#include <CGAL/Curved_kernel_3/Circle_3.h>
#include <CGAL/Circle_3.h>

#include <CGAL/Curved_kernel_3/function_objects_polynomial_sphere.h>

#include <CGAL/Curved_kernel_3/get_equation_object_on_curved_kernel_3.h>

#include <CGAL/Spherical_kernel_type_equality_wrapper.h>

namespace CGAL {
  namespace CGALi {

    template < class SphericalKernel, class LinearKernelBase >
      struct Spherical_kernel_base_ref_count: public LinearKernelBase
				 // takes classes in internal sub-namespace
      {
	typedef CGALi::Circular_arc_point_3<SphericalKernel>  Circular_arc_point_3;
        typedef CGALi::Circle_3<SphericalKernel>  Circle_3;

        // The mecanism that allows to specify reference-counting or not.
        template < typename T >
        struct Handle { typedef Handle_for<T>    type; };

        #define CGAL_Spherical_Kernel_pred(Y,Z) typedef SphericalFunctors::Y<SphericalKernel> Y; \
	    Y Z() const { return Y(); }
        #define CGAL_Spherical_Kernel_cons(Y,Z) CGAL_Spherical_Kernel_pred(Y,Z)
	
        #include <CGAL/Curved_kernel_3/interface_macros.h>
      };
    
  } // namespace CGALi

  template < class LinearKernel, class AlgebraicKernel >
    struct Spherical_kernel_3
    :  // there should be a derivation from
  // LinearKernel::Kernel_base<Self> to have types equalities for
  // the Linearkernel types
  public Spherical_kernel_type_equality_wrapper<CGALi::Spherical_kernel_base_ref_count<Spherical_kernel_3<LinearKernel,AlgebraicKernel>,
    typename LinearKernel::template Base<Spherical_kernel_3<LinearKernel,AlgebraicKernel> >::Type >,
    Spherical_kernel_3<LinearKernel,AlgebraicKernel> >
    {
      typedef Spherical_kernel_3<LinearKernel,AlgebraicKernel>      Self;

      typedef typename LinearKernel::template Base<Spherical_kernel_3<LinearKernel,AlgebraicKernel> >::Type  Linear_kernel;
      typedef AlgebraicKernel                                 Algebraic_kernel;

      //  //Please remove this if you consider it to be sloppy
      struct Curved_tag{};
      typedef Curved_tag Definition_tag;
      //  ////////////////////////////////////////////////////


      typedef typename LinearKernel::RT                       RT;
      typedef typename LinearKernel::FT                       FT;
      typedef Algebraic_kernel AK;
      typedef typename AK::Root_of_2                  Root_of_2;
      typedef typename AK::Root_for_spheres_2_3       Root_for_spheres_2_3;
      typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
      typedef typename AK::Polynomial_1_3             Polynomial_1_3;
      typedef typename AK::Polynomials_for_line_3     Polynomials_for_line_3;
      typedef std::pair< Polynomial_for_spheres_2_3, Polynomial_1_3 >
                                                 Polynomials_for_circle_3;

      // public classes
      typedef CGAL::Object Object_3;
  
    };

} // namespace CGAL

#endif // CGAL_SPHERICAL_KERNEL_3_H

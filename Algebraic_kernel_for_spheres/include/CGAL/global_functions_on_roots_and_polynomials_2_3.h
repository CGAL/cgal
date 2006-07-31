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

#ifndef CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H
#define CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H

CGAL_BEGIN_NAMESPACE

template < class AK >
inline 
typename AK::Polynomial_for_spheres_2_3
construct_polynomial_sphere_2_3(const typename AK::FT& xc,
			       const typename AK::FT& yc,
			       const typename AK::FT& zc,
			       const typename AK::FT& r_sq)
{ return AK().construct_polynomial_sphere_2_3_object()(xc, yc, zc, r_sq); }

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H

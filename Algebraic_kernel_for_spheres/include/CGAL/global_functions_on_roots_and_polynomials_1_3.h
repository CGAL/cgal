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

#ifndef CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_1_3_H
#define CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_1_3_H

CGAL_BEGIN_NAMESPACE

template < class AK >
inline 
typename AK::Polynomial_1_3
construct_polynomial_1_3(const typename AK::FT& a,
			 const typename AK::FT& b,
			 const typename AK::FT& c,
			 const typename AK::FT& d)
{ return AK().construct_polynomial_1_3_object()(a, b, c, d); }

template < class AK >
inline 
typename AK::Polynomials_for_line_3
construct_polynomials_for_line_3(const typename AK::FT& a1,
			         const typename AK::FT& b1,
			         const typename AK::FT& a2,
			         const typename AK::FT& b2,
                                 const typename AK::FT& a3,
			         const typename AK::FT& b3)
{ return AK().construct_polynomials_for_line_3_object()(a1, b1, a2, b2, a3, b3); }

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_1_3_H

// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H
#define CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H

CGAL_BEGIN_NAMESPACE

template< class AK, class OutputIterator>
inline
OutputIterator
solve( const typename AK::Polynomial_1_2 & e1,
           const typename AK::Polynomial_for_circles_2_2 & e2,
           OutputIterator res )
{ return AK().solve_object()(e1,e2,res); }

template< class AK, class OutputIterator>
inline
OutputIterator
solve( const typename AK::Polynomial_1_2 & e1,
           const typename AK::Polynomial_1_2 & e2,
           OutputIterator res )
{ return AK().solve_object()(e1,e2,res); }

template < class AK >
inline 
Sign sign_at( const typename AK::Polynomial_1_2 & equation,
	      const typename AK::Root_for_circles_2_2 r)
{ return AK().sign_at_object()(equation, r); }

template < class AK >
inline 
typename AK::Polynomial_1_2
construct_polynomial_1_2( const typename AK::RT& a,
		            const typename AK::RT& b,
		            const typename AK::RT& c)
{ return AK().construct_polynomial_1_2_object()(a, b, c); }

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H

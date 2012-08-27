// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Michael Hemmer 
// ============================================================================


#ifndef CGAL_POLYNOMIAL_FWD_H
#define CGAL_POLYNOMIAL_FWD_H

#include <CGAL/basic.h>

namespace CGAL{

template <class NT> class Polynomial; 

namespace internal{
template <class NT> inline Polynomial<NT> gcd_(const Polynomial<NT>&, const Polynomial<NT>&);
template <class NT> inline Polynomial<NT> gcd_(const Polynomial<NT>&, const Polynomial<NT>&, Field_tag);
template <class NT> inline Polynomial<NT> gcd_(const Polynomial<NT>&, const Polynomial<NT>&, Unique_factorization_domain_tag);


template <class NT> inline NT gcd_utcf_(const NT& /*a*/, const NT& /*b*/){return NT(1);}
template <class NT> inline Polynomial<NT> gcd_utcf_(const Polynomial<NT>&, const Polynomial<NT>&);
template <class NT> inline Polynomial<NT> gcd_utcf_UFD( Polynomial<NT> , Polynomial<NT>) ;
template <class NT> inline Polynomial<NT> gcd_utcf_Integral_domain(Polynomial<NT>, Polynomial<NT>);
template <class NT> inline Polynomial<NT> gcd_Euclidean_ring(Polynomial<NT>, Polynomial<NT>); 

template <class NT> inline Polynomial<NT> modular_gcd_utcf(const Polynomial<NT>&, const Polynomial<NT>&, Integral_domain_tag);
template <class NT> inline Polynomial<NT> modular_gcd_utcf(const Polynomial<NT>&, const Polynomial<NT>&, Unique_factorization_domain_tag);

// is fraction ? 
template <class NT> inline Polynomial<NT> gcd_utcf_is_fraction_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_true); 
template <class NT> inline Polynomial<NT> gcd_utcf_is_fraction_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_false);

// is type modularizable 
template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_false, Integral_domain_tag);
template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_false, Unique_factorization_domain_tag);
template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_false, Euclidean_ring_tag);

template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_true, Integral_domain_tag);
template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_true, Unique_factorization_domain_tag);
template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_true, Euclidean_ring_tag);


// template <class NT> inline NT content_utcf(const Polynomial<NT>&);
template <class NT> inline NT content_utcf_(const Polynomial<NT>&);

template <class NT, class OutputIterator1, class OutputIterator2> 
inline int filtered_square_free_factorize( Polynomial<NT>, OutputIterator1, OutputIterator2);
template <class NT,  class OutputIterator1, class OutputIterator2> 
inline int filtered_square_free_factorize_utcf( const Polynomial<NT>&, OutputIterator1, OutputIterator2);

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_utcf(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);
template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_utcf_for_regular_polynomial(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);
template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);

template <class NT> inline bool may_have_multiple_factor( const Polynomial<NT>&);
template <class NT> inline bool may_have_common_factor(const Polynomial<NT>&,const Polynomial<NT>&);

// eliminates outermost variable
template <class Coeff> 
inline Coeff resultant( 
    const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>&);
// eliminates innermost variable 
template <class Coeff> 
inline Coeff resultant_( 
    const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>&);
  


template< class Coeff >
struct Simple_matrix;

template<class NT>
internal::Simple_matrix<NT> polynomial_subresultant_matrix(
                                               CGAL::Polynomial<NT> f,
					       CGAL::Polynomial<NT> g,
					       int d=0);

template <typename Polynomial_traits_d,typename OutputIterator> inline
OutputIterator polynomial_subresultants
(typename Polynomial_traits_d::Polynomial_d A, 
 typename Polynomial_traits_d::Polynomial_d B,
 OutputIterator out);


template <typename Polynomial_traits_d,typename OutputIterator> inline
OutputIterator principal_subresultants
(typename Polynomial_traits_d::Polynomial_d A, 
 typename Polynomial_traits_d::Polynomial_d B,
 OutputIterator out);

template<typename Polynomial_traits_d,
    typename OutputIterator1, 
    typename OutputIterator2,
    typename OutputIterator3>
OutputIterator1 polynomial_subresultants_with_cofactors
(typename Polynomial_traits_d::Polynomial_d P,
 typename Polynomial_traits_d::Polynomial_d Q,
 OutputIterator1 sres_out,
 OutputIterator2 coP_out,
 OutputIterator3 coQ_out);


template <typename Polynomial_traits_d,typename OutputIterator> inline
OutputIterator
principal_sturm_habicht_sequence
(typename Polynomial_traits_d::Polynomial_d A, 
 OutputIterator out);



template<typename Polynomial_traits_d,typename OutputIterator> OutputIterator
sturm_habicht_sequence(typename Polynomial_traits_d::Polynomial_d P, 
                       OutputIterator out);

template<typename Polynomial_traits_d,
    typename OutputIterator1,
    typename OutputIterator2,
    typename OutputIterator3> 
OutputIterator1
sturm_habicht_sequence_with_cofactors
(typename Polynomial_traits_d::Polynomial_d P,
       OutputIterator1 out_stha,
 OutputIterator2 out_f,
 OutputIterator3 out_fx);

} // namespace internal


} // namespace CGAL


#include <CGAL/Polynomial.h>

#endif

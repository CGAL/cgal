// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//
// ============================================================================


#ifndef CGAL_POLYNOMIAL_FWD_H
#define CGAL_POLYNOMIAL_FWD_H

#include <CGAL/basic.h>

namespace CGAL{

template <class NT> class Polynomial; 

namespace CGALi{
template <class NT> inline Polynomial<NT> gcd (const Polynomial<NT>&, const Polynomial<NT>&);
template <class NT> inline Polynomial<NT> gcd_(const Polynomial<NT>&, const Polynomial<NT>&, Field_tag);
template <class NT> inline Polynomial<NT> gcd_(const Polynomial<NT>&, const Polynomial<NT>&, Unique_factorization_domain_tag);


template <class NT> inline NT gcd_utcf_(const NT& a, const NT& b){return NT(1);}
template <class NT> inline Polynomial<NT> gcd_utcf(const Polynomial<NT>&, const Polynomial<NT>&);

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
template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_false, Field_tag);

template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_true, Integral_domain_tag);
template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_true, Unique_factorization_domain_tag);
template <class NT> inline Polynomial<NT> gcd_utcf_modularizable_algebra_( const Polynomial<NT>&, const Polynomial<NT>&, ::CGAL::Tag_true, Field_tag);


// template <class NT> inline NT content_utcf(const Polynomial<NT>&);
template <class NT> inline NT content_utcf_(const Polynomial<NT>&);

template <class NT, class OutputIterator1, class OutputIterator2> 
inline int filtered_square_free_factorization( Polynomial<NT>, OutputIterator1, OutputIterator2);
template <class NT,  class OutputIterator1, class OutputIterator2> 
inline int filtered_square_free_factorization_utcf( const Polynomial<NT>&, OutputIterator1, OutputIterator2);

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorization_utcf(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);
template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorization_utcf_for_regular_polynomial(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorization(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);
template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorization_for_regular_polynomial(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);

template <class NT> inline bool may_have_multiple_factor( const Polynomial<NT>&);
template <class NT> inline bool may_have_common_factor(const Polynomial<NT>&,const Polynomial<NT>&);


} // namespace CGALi
} // namespace CGAL


#include <CGAL/Polynomial.h>

#endif

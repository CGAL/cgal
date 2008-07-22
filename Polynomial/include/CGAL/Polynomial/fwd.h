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

template< class Coeff >
struct Simple_matrix;

template<class NT>
CGALi::Simple_matrix<NT> polynomial_subresultant_matrix(
                                               CGAL::Polynomial<NT> f,
					       CGAL::Polynomial<NT> g,
					       int d=0);


template <typename OutputIterator, typename NT> inline
    OutputIterator polynomial_subresultants(CGAL::Polynomial<NT> A, CGAL::Polynomial<NT> B,
                                            OutputIterator out);
template <typename OutputIterator, typename NT> inline
    OutputIterator principal_subresultants(CGAL::Polynomial<NT> A, CGAL::Polynomial<NT> B,
                                            OutputIterator out);

template <typename OutputIterator, typename NT > inline
      OutputIterator principal_sturm_habicht_sequence(CGAL::Polynomial<NT> A, 
                                                      OutputIterator out);

template<typename NT,
    typename OutputIterator1, 
    typename OutputIterator2,
    typename OutputIterator3>
    OutputIterator1 prs_subresultants_with_cofactors(CGAL::Polynomial<NT> P,
                                                     CGAL::Polynomial<NT> Q,
                                                     OutputIterator1 sres_out,
                                                     OutputIterator2 coP_out,
                                                     OutputIterator3 coQ_out);


template<typename OutputIterator, typename NT> OutputIterator
    sturm_habicht_sequence(CGAL::Polynomial<NT> P, OutputIterator out);

  template<typename OutputIterator1,
    typename OutputIterator2,
    typename OutputIterator3,
    typename NT> OutputIterator1
    sturm_habicht_sequence_with_cofactors(CGAL::Polynomial<NT> P, 
                                          OutputIterator1 out_stha,
                                          OutputIterator2 out_f,
                                          OutputIterator3 out_fx);


} // namespace CGALi



} // namespace CGAL


#include <CGAL/Polynomial.h>

#endif

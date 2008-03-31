// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_POLYNOMIAL_RESULTANT_H
#define CGAL_POLYNOMIAL_RESULTANT_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial/prs_resultant.h>
#include <CGAL/Polynomial/bezout_matrix.h>

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

    template <class NT> inline
    NT resultant_(Polynomial<NT> A, Polynomial<NT> B) {
        typedef typename Algebraic_structure_traits<NT>::Algebraic_category Algebraic_category;
        return resultant_(A,B,Algebraic_category());     
    }
    
    // the general function for CGAL::Integral_domain_without_div_tag
    template < class NT> inline 
    NT resultant_(Polynomial<NT> A, 
                  Polynomial<NT> B,
                  Integral_domain_without_division_tag){
        return hybrid_bezout_subresultant(A,B,0);
    }
    
    // the specialization for CGAL::Unique_factorization_domain_tag
    template < class NT> inline
    NT resultant_(Polynomial<NT> A, 
                  Polynomial<NT> B,
                  Unique_factorization_domain_tag ){
        return prs_resultant_ufd(A,B);
    }
        
    template <class NT> inline
    NT resultant_field(Polynomial<NT> A, Polynomial<NT> B, ::CGAL::Tag_false) {       
        return prs_resultant_field(A,B);       
    } 
    template <class NT> inline
    NT resultant_field(Polynomial<NT> A, Polynomial<NT> B, ::CGAL::Tag_true) {       
        return prs_resultant_decompose(A,B);       
    } 
    
    // the specialization for CGAL::Field_tag
    template < class NT> inline
    NT resultant_(Polynomial<NT> A, Polynomial<NT> B, Field_tag){
        // prs_resultant_decompose can only be used if 
        // NT is decomposable && the resulting type is at least a UFDomain
        typedef typename Fraction_traits<NT>::Is_fraction Is_decomposable;
        const bool is_decomposable 
            = ::boost::is_same<Is_decomposable, ::CGAL::Tag_true>::value;
        
        
        typedef typename Fraction_traits<NT>::Numerator_type NUM;
        typedef typename Algebraic_structure_traits<NUM>::Algebraic_category ALG_TYPE;
        const bool is_inheritance 
            = ::boost::is_base_and_derived<Unique_factorization_domain_tag,ALG_TYPE>::value ||
              ::boost::is_same<Unique_factorization_domain_tag,ALG_TYPE>::value;

        return resultant_field(
          A, B,
          typename ::boost::mpl::if_c< is_decomposable && is_inheritance,
                                 ::CGAL::Tag_true,
                                 ::CGAL::Tag_false >::type() );
    }    

    /*! \ingroup CGAL_Polynomial
     *  \relates CGAL::Polynomial
     *  \brief compute the resultant of the polynomials \c A and \c B
     *
     *  The way the resultant is computed depends on the Algebraic_category. 
     *  In general the resultant will be computed by the function
     *  CGAL::hybrid_bezout_subresultant, but if possible the function
     *  CGAL::prs_resultant_ufd or CGAL::prs_resultant_field are used. 
     *  
     *  Up to now it is not clear, that the functions based on the polynomial
     *  remainder sequence are faster than the one based on the bezoutian. 
     *  Thus you can use CGAL::hybrid_bezout_subresultant instead, which will 
     *  work for any Algebraic_category
     */
    template <class NT> inline
    NT resultant(Polynomial<NT> A, Polynomial<NT> B) {
        return CGALi::resultant_(A, B);
    }   

} // namespace CGALi

CGAL_END_NAMESPACE

#endif// CGAL_POLYNOMIAL_RESULTANT_H

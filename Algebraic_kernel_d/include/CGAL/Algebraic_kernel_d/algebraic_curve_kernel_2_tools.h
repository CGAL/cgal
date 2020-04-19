// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================


#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_2_TOOLS
#define CGAL_ALGEBRAIC_CURVE_KERNEL_2_TOOLS 1

#include <iterator>

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_d/enums.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Float_traits.h>
#include <CGAL/convert_to_bfi.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <boost/numeric/interval.hpp>
#include <CGAL/Algebraic_kernel_d/bound_between_1.h>
#include <CGAL/Coercion_traits.h>

namespace CGAL {

namespace internal {

/*
 * \brief Function for merging two sets
 *
 * This function is similar to the \c std::union_set operation.
 * Additionally, it provides a sequence of CGAL::internal::Three_valued
 * providing information to which input set the corresponding root
 * in the merged sequence belonged
 *
 * The BinaryFunction must have the Result type CGAL::Comparison_result.
 */
template<typename BinaryFunction,
      typename InputIterator1,typename InputIterator2,
      typename OutputIterator1,typename OutputIterator2>
std::pair<OutputIterator1,OutputIterator2>
set_union_with_source(InputIterator1 first_begin,
                      InputIterator1 first_end,
                      InputIterator2 second_begin,
                      InputIterator2 second_end,
                      OutputIterator1 merged_values,
                      OutputIterator2 merged_values_info,
                      BinaryFunction compare) {

    InputIterator1 first_it=first_begin;
    InputIterator2 second_it=second_begin;

    while((first_it != first_end) || (second_it!=second_end)) {
        if(first_it == first_end) {
            *merged_values=*second_it++;
            ++merged_values;
            *merged_values_info++ = CGAL::internal::ROOT_OF_SECOND_SET;
            continue;
        }
        if(second_it == second_end) {
            *merged_values++=*first_it++;
            *merged_values_info++ = CGAL::internal::ROOT_OF_FIRST_SET;
            continue;
        }

        CGAL::Comparison_result c = compare(*first_it,*second_it);

        if(c==CGAL::EQUAL) {
            *merged_values++=*first_it++;
            ++second_it;
            *merged_values_info++ = CGAL::internal::ROOT_OF_BOTH_SETS;
            continue;
        }
        if(c==CGAL::SMALLER) {
            *merged_values++=*first_it++;
            *merged_values_info++ = CGAL::internal::ROOT_OF_FIRST_SET;
            continue;
        }
        if(c==CGAL::LARGER) {
            *merged_values++=*second_it++;
            *merged_values_info++ = CGAL::internal::ROOT_OF_SECOND_SET;
            continue;
        }
    }
    return std::make_pair(merged_values,merged_values_info);
}

/*
 * \brief Removes the leading term of the polynomial \c f as long as it
 * vanishes at \c alpha
 *
 */
template<typename Algebraic_kernel_d_1,typename Poly_2, typename Algebraic_real>
Poly_2 poly_non_vanish_leading_term(Algebraic_kernel_d_1* kernel,
                                    const Poly_2& pol,
                                    Algebraic_real alpha) {
    Poly_2 f(pol);
    while(true) {
        if(kernel->is_zero_at_1_object()
           (CGAL::leading_coefficient(f),alpha)) {
            typename Poly_2::const_iterator poly_end = f.end();
            if(f.begin()==poly_end) {
                break;
            }
            poly_end--;
            f=Poly_2(f.begin(),poly_end);
        }
        else {
            break;
        }
    }
    return f;
}

/*!
 * \brief finds a Rational value left of an Algebraic real alpha
 */
template<typename AlgebraicKernel_1> typename AlgebraicKernel_1::Bound
  bound_left_of(const AlgebraicKernel_1* kernel,
                typename AlgebraicKernel_1::Algebraic_real_1 ar) {

    typedef AlgebraicKernel_1 Algebraic_kernel_d_1;

    typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;
    typedef typename Algebraic_kernel_d_1::Bound Bound;

    switch( CGAL::sign( ar ) ) {
    case(CGAL::ZERO): {
        return Bound(-1);
        break;
    }
    case(CGAL::POSITIVE): {
        return Bound(0);
        break;
    }
    case(CGAL::NEGATIVE): {
        Algebraic_real_1 small_value
          = kernel->construct_algebraic_real_1_object()
          (Bound(2)*kernel->approximate_absolute_1_object()(ar,1).first);

        return kernel->bound_between_1_object()(small_value,ar);
        // = small_value.rational_between(ar);
        //= ar.low()-1;
    }
    }
    // never reached
    return Bound(0);
}

/*!
 * \brief finds a Rational value rightt of an Algebraic real alpha
 */
template<typename AlgebraicKernel_1> typename AlgebraicKernel_1::Bound
  bound_right_of(const AlgebraicKernel_1* kernel,
                 typename AlgebraicKernel_1::Algebraic_real_1 ar) {

    return -bound_left_of(kernel,-ar);

}


/*!
 * \brief Produces intermediate rational values for a list of
 * algebraic reals.
 *
 * For a list of Algebraic real values with \c n elements, a list with
 * <tt>n+1</tt> elements of rational values is given such that the
 * <tt>i</tt>th element is
 * between the <tt>i</tt>th and the <tt>(i+1)</tt>th element of the input list
 *
 * The input list must be in increasing order
 */
template<typename AlgebraicKernel_1,
         typename InputIterator,
         typename OutputIterator>
  OutputIterator find_intermediate_values(const AlgebraicKernel_1* kernel,
                                          InputIterator start,
                                          InputIterator end,
                                          OutputIterator output) {
    CGAL_static_assertion
      ((::boost::is_same
        <typename AlgebraicKernel_1::Algebraic_real_1,
        typename std::iterator_traits<InputIterator>::value_type >::value));

    typedef typename AlgebraicKernel_1::Bound Bound;
    if(start==end) {
        // Empty vector, create one element
        *output++=Bound(0);
        return output;
    }
    *output++=bound_left_of(kernel,*start);

    InputIterator it_1(start),it_2(start);
    ++it_2;
    while(it_2 != end) {
        CGAL_assertion(it_1->compare(*it_2)==CGAL::SMALLER);
        Bound beta
          = kernel->bound_between_1_object()(*it_1,*it_2);
        *output++=beta;
        ++it_1;
        ++it_2;
    }
    *output++=bound_right_of(kernel,*it_1);

    return output;
}



// Used internally for zero_test_bivariate

namespace for_zero_test_bivariate {

template<typename Poly_coer_1,typename Polynomial_1>
  void cast_back_utcf(const Poly_coer_1& p,Polynomial_1& q) {
  // We can assume that both template arguments are polynomial types
  typedef CGAL::Fraction_traits<Poly_coer_1> FT;
  CGAL_static_assertion((::boost::is_same<typename FT::Is_fraction,
                       CGAL::Tag_true>::value));
  typedef typename FT::Numerator_type Numerator;
  typedef typename FT::Denominator_type Denominator;
  typedef CGAL::Coercion_traits<Numerator,Polynomial_1> Num_coercion;
  CGAL_static_assertion((::boost::is_same
                       <Polynomial_1,
                       typename Num_coercion::Type>::value));
  Numerator p_num;
  Denominator p_denom;
  typename FT::Decompose()(p,p_num,p_denom);
  q = typename Num_coercion::Cast()(p_num);
}

template<typename A> void cast_back_utcf(const A& p, A& q) {
  q = p;
}

} // of namespace for_zero_test_bivariate

/*
 * \brief  Symbolic zero test.
 *
 * Checks whether <tt>h(x,y(x))=0</tt>, where <tt>y(x)</tt> is a rational
 * expression in terms of \c x, i.e. <tt>y=p/q</tt> with <tt>p,q</tt>
 * univariate polynomials
 */
template<typename AlgebraicCurveKernel_2>
  bool zero_test_bivariate
  (const AlgebraicCurveKernel_2* kernel,
   const typename AlgebraicCurveKernel_2::Algebraic_real_1& alpha,
   const typename AlgebraicCurveKernel_2::Polynomial_2& h,
   const typename AlgebraicCurveKernel_2::Polynomial_1& p,
   const typename AlgebraicCurveKernel_2::Polynomial_1& q) {

    bool result;
    typedef typename AlgebraicCurveKernel_2::Polynomial_1 Polynomial_1;
#if !CGAL_ACK_USE_NO_REDUCTION_MODULO_RESULTANT

    //typedef typename AlgebraicCurveKernel_2::Algebraic_real_1 Algebraic_real_1;
    typedef typename AlgebraicCurveKernel_2::Bound Bound;
    typedef typename AlgebraicCurveKernel_2::Coefficient Coefficient;
    typedef typename AlgebraicCurveKernel_2::Polynomial_2 Polynomial_2;

    typedef CGAL::Coercion_traits<Bound,Coefficient> Coercion;
    typedef typename Coercion::Type Coercion_type;
    typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
      ::template Rebind<Coercion_type,1>::Other::Type Poly_coer_1;

    typename Coercion::Cast cast;

    bool general = ! alpha.is_rational();

    Poly_coer_1 p_rat = typename CGAL::Polynomial_traits_d<Poly_coer_1>
      ::Construct_polynomial()
      (boost::make_transform_iterator
       (p.begin(),cast),
       boost::make_transform_iterator
       (p.end(),cast));
    Poly_coer_1 q_rat = typename CGAL::Polynomial_traits_d<Poly_coer_1>
      ::Construct_polynomial()
      (boost::make_transform_iterator
       (q.begin(),cast),
       boost::make_transform_iterator
       (q.end(),cast));

    if(general) {


        Poly_coer_1 modulus = typename CGAL::Polynomial_traits_d<Poly_coer_1>
        ::Construct_polynomial()
            (boost::make_transform_iterator
             (alpha.polynomial().begin(),cast),
             boost::make_transform_iterator
             (alpha.polynomial().end(),cast));

/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << "Mod: " << modulus << std::endl;
  #endif
*/
        p_rat=CGAL::mod(p_rat,modulus);
        q_rat=CGAL::mod(q_rat,modulus);

        int n = CGAL::degree(h,1);
        // Create the powers of p and q mod modulus
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "precomp powers.." << std::flush;
#endif
*/
        std::vector<Poly_coer_1> p_powers(n+1),q_powers(n+1);
        p_powers[0]=Poly_coer_1(Bound(1));
        q_powers[0]=Poly_coer_1(Bound(1));
        Poly_coer_1 intermediate;
        for(int i=1;i<=n;i++) {
/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << i << ": mult.." << std::flush;
  #endif
*/
            intermediate=p_powers[i-1]*p_rat;
/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << "mod.." << std::flush;
  #endif
*/
            p_powers[i]=CGAL::mod(intermediate,modulus);
/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << "simpl.." << std::flush;
  #endif
*/
            p_powers[i].simplify_coefficients();
/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << "mult.." << std::flush;
  #endif
*/
            intermediate=q_powers[i-1]*q_rat;
/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << "mod.." << std::flush;
  #endif
*/
            q_powers[i]=CGAL::mod(intermediate,modulus);
/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << "simpl.." << std::flush;
  #endif
*/
            q_powers[i].simplify_coefficients();
        }
/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << "done\ncomp rat pol.." << std::flush;
  #endif
*/

        Poly_coer_1 curr_coeff,curr_fac;
        Poly_coer_1 h_0_rat(Coercion_type(0));
        for(int i=0;i<=n;i++) {
          Poly_coer_1 tmp_pol = typename CGAL::Polynomial_traits_d<Poly_coer_1>
            ::Construct_polynomial()
                (boost::make_transform_iterator
                 (h[i].begin(),cast),
                 boost::make_transform_iterator
                 (h[i].end(),cast));
            curr_fac=CGAL::mod
                (tmp_pol*p_powers[i]*q_powers[n-i], modulus);
            h_0_rat+=curr_fac;
        }

        Polynomial_1 h_0_utcf;
        for_zero_test_bivariate::cast_back_utcf(h_0_rat,h_0_utcf);

        return kernel->is_zero_at_1_object() (h_0_utcf,alpha);
    }
    else {
        Coercion_type b = cast(alpha.rational()),
          p_b=CGAL::evaluate(p_rat,b),q_b=CGAL::evaluate(q_rat,b);
        int n = CGAL::degree(h,1);
        Coercion_type eval(0);
        for(int i=0;i<=n;i++) {
          Poly_coer_1 h_i_rat = typename CGAL::Polynomial_traits_d<Poly_coer_1>
            ::Construct_polynomial()
            (boost::make_transform_iterator
             (CGAL::get_coefficient(h,i).begin(),cast),
             boost::make_transform_iterator
             (CGAL::get_coefficient(h,i).end(),cast));
          eval+=CGAL::evaluate(h_i_rat,b)
            *CGAL::ipower(p_b,i)*CGAL::ipower(q_b,n-i);
        }
        result=(CGAL::sign(eval)==CGAL::ZERO);
    }
#else
#warning Uses no reduction modulo resultant!
    Polynomial_1 h_0=CGAL::evaluate_homogeneous(h,p,q);
    result= kernel->is_zero_at_1_object() (h_0,alpha);
#endif

    return result;

}



} // namespace internal


} //namespace CGAL

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_2_TOOLS

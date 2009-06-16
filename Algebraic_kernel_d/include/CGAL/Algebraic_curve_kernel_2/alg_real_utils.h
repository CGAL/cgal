// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================


#ifndef CGAL_ALG_REAL_UTILS
#define CGAL_ALG_REAL_UTILS 1

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

#include <CGAL/Coercion_traits.h>

// Workaround to solve namepsace-problems for boost when applying the minimum
// of Sqrt_extension

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/* 
 * \brief Function for merging two sets
 *
 * This function is similar to the \c std::union_set operation.
 * Additionally, it provides a sequence of AcX::Three_valued
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
            *merged_values_info++ = CGAL::CGALi::ROOT_OF_SECOND_SET;
            continue;
	}
	if(second_it == second_end) {
            *merged_values++=*first_it++;
            *merged_values_info++ = CGAL::CGALi::ROOT_OF_FIRST_SET;
            continue;
	}

	CGAL::Comparison_result c = compare(*first_it,*second_it);

	if(c==CGAL::EQUAL) {
            *merged_values++=*first_it++;
            ++second_it;
            *merged_values_info++ = CGAL::CGALi::ROOT_OF_BOTH_SETS;
            continue;
	}
	if(c==CGAL::SMALLER) {
            *merged_values++=*first_it++;
            *merged_values_info++ = CGAL::CGALi::ROOT_OF_FIRST_SET;
            continue;
	}
	if(c==CGAL::LARGER) {
            *merged_values++=*second_it++;
            *merged_values_info++ = CGAL::CGALi::ROOT_OF_SECOND_SET;
            continue;
	}
    }
    return std::make_pair(merged_values,merged_values_info);
}



/*! \brief Tries to find a SIMPLE rational q with a<q<b.
 *
 * In this context, simple means that the denominator of <tt>q</tt>
 * is a power of two, and is not too big. There is no guarantee to find
 * the rational value between <tt>a</tt> and <tt>b</tt> of minimal
 * bit size.
 */
template<typename Algebraic_real>
typename Algebraic_real::Rational
simple_rational_between(const Algebraic_real& a,
                        const Algebraic_real&b) {

    //srb.start();
    typedef typename Algebraic_real::Rational Rational;
    typename CGAL::Fraction_traits<Rational>::Compose compose;
    typedef typename 
        CGAL::Get_arithmetic_kernel<Rational>::Arithmetic_kernel AK;
    typedef typename AK::Bigfloat_interval Bigfloat_interval;
    typedef typename CGAL::Bigfloat_interval_traits<Bigfloat_interval>
        ::Boundary Bigfloat;
    typedef typename AK::Integer Integer;

    long old_prec = CGAL::get_precision(Bigfloat_interval());

    CGAL_assertion(a<b);
    //std::cout << "Intermediate1: " << CGAL::to_double(a) << " " << CGAL::to_double(b) << std::endl;
    /*
     * First, refine a and b until their isolating intervals are disjoint
     * Therefore, the bigger interval is refined in each substep
     */
    //srb_a.start();
    if(a.high() >= b.low()) {
        Rational size_a=a.high()-a.low(),
            size_b=b.high() - b.low();
        while(a.high() >= b.low()) {
            if(size_a < size_b) {
                b.refine();
                size_b=b.high() - b.low();
            } else {
                a.refine();
                size_a=a.high()-a.low();
            }
        }
    }
    //srb_a.stop();

    //srb_b.start();
    Bigfloat x=CGAL::upper(CGAL::convert_to_bfi(a.high())),
        y=CGAL::lower(CGAL::convert_to_bfi(b.low()));
    if(x>=y) {
        Rational size_a=a.high() - a.low(),
            size_b=b.high() - b.low(),
            size_max = size_a>size_b ? size_a : size_b,
            size_int = b.low()-a.high();
        while(x>=y) {
            //std::cout << "x and y: " <<  x << " and " << y << std::endl;
            //std::cout << "sizes: " << CGAL::to_double(size_int) << " " << CGAL::to_double(size_max) << std::endl;
            if(size_int>size_max) {
                CGAL::set_precision(Bigfloat_interval(),
                                   2*CGAL::get_precision(Bigfloat_interval()));
                x=CGAL::upper(CGAL::convert_to_bfi(a.high()));
                y=CGAL::lower(CGAL::convert_to_bfi(b.low()));
            } else {
                if(size_a < size_b) {
                    b.refine();
                    size_b=b.high() - b.low();
                    y=CGAL::lower(CGAL::convert_to_bfi(b.low()));
                } else {
                    a.refine();
                    size_a=a.high()-a.low();    
                    x=CGAL::upper(CGAL::convert_to_bfi(a.high()));
                }
                size_max = size_a>size_b ? size_a : size_b;
                size_int = b.low()-a.high();
            }
        }
    }
    //srb_b.stop();
    //std::cout << "Intermediate2: " << x << " " << y << std::endl;
    typename CGAL::CGALi::Float_traits<Bigfloat>::Get_mantissa mantissa;
    typename CGAL::CGALi::Float_traits<Bigfloat>::Get_exponent exponent;

    Integer x_m = mantissa(x),
        y_m=mantissa(y);
    long x_e = exponent(x),  y_e = exponent(y);
    //std::cout << "Floats1: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;
    if((x_m > 0 &&  y_m < 0) || x_m < 0 && y_m > 0) {
        //srb.stop();
        return Rational(0);
    }
    bool negative=false;
    if(x_m<=0 && y_m <=0) {
        x_m=-x_m;
        y_m=-y_m;
        std::swap(x_m,y_m);
        std::swap(x_e,y_e);
        negative=true;
    }
    // Now, we have that (x_m,x_e) represents a number smaller than (y_m,y_e)
    //srb_c.start();
    //std::cout << "Floats2: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;

    // As long as the mantissa is even, simplify
    while(x_m != 0 && (x_m & 1)==0 ) {
        x_m=x_m >> 1;
        x_e++;
    }
    while(y_m != 0 && (y_m & 1)==0 ) {
        y_m=y_m >> 1;
        y_e++;
    }
    //srb_c.stop();
    //std::cout << "Floats3: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;

    // Bring both numbers to a common exponent
    //srb_d.start();
    long min_e = x_e < y_e ? x_e : y_e;
    while(x_e > min_e) {
        x_m=x_m << 1;
        x_e--;
    }
    while(y_e > min_e) {
        y_m=y_m << 1;
        y_e--;
    }
    //srb_d.stop();
    CGAL_assertion(y_e==x_e && x_e==min_e);

    //std::cout << "Floats4: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;

    // Avoid mantissas to have difference one
    if(y_m-x_m==Integer(1)) {
        x_m=x_m << 1;
        y_m=y_m << 1;
        x_e--;
        y_e--;
        min_e--;
    }
    //std::cout << "Floats5: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;
    Integer final_mantissa(0);
    //srb_e.start();
    long x_log = x_m==Integer(0) ? -1 : CGAL::CGALi::floor_log2_abs(x_m),
        y_log = y_m==Integer(0) ? -1 : CGAL::CGALi::floor_log2_abs(y_m),
        old_log = y_log;
    while(x_log==y_log) {
        //std::cout << "here" << std::endl;
        while(old_log > y_log) {
            final_mantissa = final_mantissa << 1;
            old_log--;
        }
        CGAL_assertion((x_m & ((Integer(1) << x_log) - 1)) == x_m - CGAL::ipower(Integer(2),x_log));
        x_m = x_m & ((Integer(1) << x_log) - 1); // x_m - CGAL::ipower(Integer(2),x_log);
        y_m = y_m & ((Integer(1) << y_log) - 1); // y_m - CGAL::ipower(Integer(2),y_log);

        final_mantissa++;
        old_log=y_log;
        x_log = x_m==0 ? -1 : CGAL::CGALi::floor_log2_abs(x_m);
        y_log = y_m==0 ? -1 : CGAL::CGALi::floor_log2_abs(y_m);
    }
    //srb_e.stop();
    // Now, x_log != y_log, in fact, y_log is greater
    CGAL_assertion(x_log<y_log);
    //srb_f.start();
    while(old_log > y_log) {
        final_mantissa = final_mantissa << 1;
        old_log--;
    }
    if((y_m & ((Integer(1) << y_log) - 1 ))==0) { // y_m - CGAL::ipower(Integer(2),y_log)==0) {
        // Now, the constructed value would be equal to
        while(y_log!=0 && x_log==y_log-1) {
            final_mantissa = final_mantissa << 1;
            final_mantissa++;
            y_log--;
            x_m = x_m==0 ? 0 : x_m & ((Integer(1) << x_log) - 1); //x_m - CGAL::ipower(Integer(2),x_log);
            x_log = x_m==0 ? -1 : CGAL::CGALi::floor_log2_abs(x_m);
        }
        final_mantissa = final_mantissa << 1;
        final_mantissa++;
        y_log--;
    } else {
        final_mantissa++;
    }
    //srb_f.stop();
    min_e += y_log;
    Rational rat_between;
    //std::cout << "Min_e: " << min_e << std::endl;
    if(min_e > 0) {
        rat_between = compose(final_mantissa << min_e,
                              Integer(1));
    } else {
        rat_between = compose(final_mantissa, Integer(1) << -min_e);
    }
    if(negative) {
        rat_between = -rat_between;
    }
    //std::cout << "Result: " << a.high() << " " << rat_between << " " << b.low() << std::endl;
    CGAL_assertion(a.high() < rat_between);
    CGAL_assertion(b.low() > rat_between);
    CGAL::set_precision(Bigfloat_interval(),old_prec);
    //srb.stop();
    return rat_between;
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
template<typename InputIterator, typename OutputIterator> 
OutputIterator find_intermediate_values(InputIterator start,
                                        InputIterator end,
                                        OutputIterator output) {
    typedef typename std::iterator_traits<InputIterator>::value_type 
	Alg_real;
    typedef typename Alg_real::Rational Rational;
    if(start==end) {
	// Empty vector, create one element
	*output=Rational(0);
	++output;
	return output;
    }
    switch( CGAL::sign( start->low() ) ) {
    case(CGAL::ZERO): {
        *output = Rational(-1);
        break;
    }
    case(CGAL::POSITIVE): {
        *output = Rational(0);
        break;
    }
    case(CGAL::NEGATIVE): {
        Alg_real small_value = Alg_real(Rational(2)*start->low());
        *output 
            = simple_rational_between(small_value,*start);
    }
    }
    output++;
      
    ++output;
    InputIterator it_1(start),it_2(start);
    ++it_2;
    while(it_2 != end) {
	CGAL_assertion(it_1->compare(*it_2)==CGAL::SMALLER);
	Rational beta 
            = simple_rational_between(*it_1,*it_2);
	*output=beta;
	++output;
	++it_1;
	++it_2;
    }
    switch( CGAL::sign( it_1->high() ) ) {
    case(CGAL::ZERO): {
        *output = Rational(1);
        break;
    }
    case(CGAL::NEGATIVE): {
        *output = Rational(0);
        break;
    }
    case(CGAL::POSITIVE): {
        Alg_real big_value = Alg_real(Rational(2)*it_1->high());
        *output 
            = simple_rational_between(*it_1,big_value);
    }
    }
    output++;
    return output;
}

/*!
 * \brief finds a Rational value left of an Algebraic real alpha
 */
template<typename AlgebraicReal> typename AlgebraicReal::Rational 
simple_rational_left_of(AlgebraicReal ar) {

    typedef AlgebraicReal Algebraic_real;
    typedef typename Algebraic_real::Rational Rational;

    switch( CGAL::sign( ar ) ) {
    case(CGAL::ZERO): {
        return Rational(-1);
        break;
    }
    case(CGAL::POSITIVE): {
        return Rational(0);
        break;
    }
    case(CGAL::NEGATIVE): {
        Algebraic_real small_value 
            = Algebraic_real(Rational(2)*ar.low());
    
        return simple_rational_between(small_value,ar);
        // = small_value.rational_between(ar);
        //= ar.low()-1;
    }
    }
    // never reached
    return Rational(0);
}

/*!
 * \brief finds a Rational value rightt of an Algebraic real alpha
 */
template<typename AlgebraicReal> typename AlgebraicReal::Rational 
simple_rational_right_of(AlgebraicReal ar) {
    
    typedef AlgebraicReal Algebraic_real;
    typedef typename Algebraic_real::Rational Rational;

    return -simple_rational_left_of(-ar);

}
  
/*!
 * \brief Tries to get the (non-zero) sign of <tt>f(alpha)</tt> with interval
 * arithmetic.
 *
 * Performs refinements of alpha and evaluates \c f at \c alpha. If
 * zero is still in the interval when the isolating interval of \c alpha
 * is smaller then \c eps, CGAL::ZERO
 * is returned, denoting that the sign could not be computed. Otherwise,
 * the correct sign of \c f(alpha) is returned
 *
 * If \c n is set to zero (default), \c alpha is refined 
 * until the sign can be determined correctly. This may lead to an 
 * infinite loop, if <tt>f(alpha)=0</tt>
 *
 */
template<typename NT,typename Polynomial_1>
CGAL::Sign estimate_sign_at(NT alpha,
                            const Polynomial_1& f,
                            long max_precision=0) {
    typedef CGAL::CGALi::Bitstream_coefficient_kernel<NT> 
        Bitstream_coefficient_kernel;
    typedef typename Bitstream_coefficient_kernel::Bigfloat_interval BFI;
    long old_prec = CGAL::get_precision(BFI());
    long prec = 16;

    typename Bitstream_coefficient_kernel::Convert_to_bfi convert_to_bfi;

    CGAL::Sign sign=CGAL::ZERO;

    while(max_precision==0 || prec<=max_precision) {
        CGAL::set_precision(BFI(),prec);
        BFI eval = f.evaluate(convert_to_bfi(alpha));
        if(! CGAL::in_zero(eval)) {
            sign=CGAL::sign(eval);
            CGAL_assertion(sign!=CGAL::ZERO);
            break;
        } else {
            prec*=2;
        }
    }
    CGAL::set_precision(BFI(),old_prec);
    return sign;
}
     

/*!
 * \brief Evaluates poynomial at intervals with the Horner scheme
 */
template<typename Poly1_,typename Boundary_>
boost::numeric::interval
<typename CGAL::Coercion_traits<typename Poly1_::NT, Boundary_>::Type >
evaluate_iv(Poly1_ f,boost::numeric::interval<Boundary_> iv) {

    typedef Boundary_ Boundary;
    typedef Poly1_ Poly1;
    typedef typename 
        CGAL::Coercion_traits<typename Poly1::NT, Boundary_>::Type Coercion;
	//typename CGAL::Coercion_traits<typename Poly1::NT, Boundary_>::Cast cast;
    CGAL_assertion(f.degree()>=0);
    typename CGAL::Coercion_traits<typename Poly1::NT, Boundary_>::Cast cast;
    typedef boost::numeric::interval<Coercion> Coercion_interval;
    Coercion_interval iv_cast(cast(iv.lower()),cast(iv.upper()));
    int n=f.degree();
    Coercion_interval ret(cast(f[n]),cast(f[n]));
    for(int i=n-1;i>=0;i--) {
        ret *= iv_cast;
        ret += Coercion_interval(cast(f[i]),cast(f[i]));
    }
    return ret;
}

//! \brief replaces the \c NiX::is_root_of function, using \c AcX::ntl_gcd
template<typename Polynomial_2,typename AlgebraicReal>
bool is_root_of(const AlgebraicReal & x,Polynomial_2 p) {
    if(x.is_rational()) {
        typedef typename AlgebraicReal::Rational Rational;
        Rational exact_value = x.rational();
        typedef typename CGAL::Coercion_traits<Rational, 
            typename Polynomial_2::NT>::Type
            Type;
        return CGAL::sign(p.evaluate(exact_value))==CGAL::ZERO; 
    }
    else {
        CGAL_precondition(p.degree()>=0);
      
        if(p.degree()==0) return p.is_zero();
        Polynomial_2 g=typename CGAL::Polynomial_traits_d<Polynomial_2>
            ::Gcd_up_to_constant_factor()(p,x.polynomial());
        return g.sign_at(x.low())!=g.sign_at(x.high());
      
    }
    return false; 
}

/*
 * \brief Removes the leading term of the polynomial \c f as long as it
 * vanishes at \c alpha
 *
 */
template<typename Poly_2, typename Algebraic_real>
Poly_2 poly_non_vanish_leading_term(const Poly_2& pol,Algebraic_real alpha) {
    Poly_2 f(pol);
    while(true) {
	if(CGAL::CGALi::is_root_of(alpha,f.lcoeff())) {
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
	  
} // namespace CGALi


CGAL_END_NAMESPACE

#endif // CGAL_ALG_REAL_UTILS

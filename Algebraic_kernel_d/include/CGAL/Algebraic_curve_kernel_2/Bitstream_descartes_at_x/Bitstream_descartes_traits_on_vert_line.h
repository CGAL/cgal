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

#ifndef CGAL_BITSTREAM_DESCARTES_TRAITS_ON_VERT_LINE
#define CGAL_BITSTREAM_DESCARTES_TRAITS_ON_VERT_LINE 1

#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/alg_real_utils.h>
#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/Best_approximation_cache.h>

CGAL_BEGIN_NAMESPACE 

    // Workaround for min/max. This is just temporary!
    // Necessary to prevent compiler confusion about std::min/max 
    // and CGAL::min/max. The better solution is to put the CGAL functions into
    // some seperate namespace and include it into CGAL-namepsace with the
    // using command. The workaround works only for Sqrt_extensions
template< class NT, class Root >
inline Sqrt_extension<NT,Root> min BOOST_PREVENT_MACRO_SUBSTITUTION
( const Sqrt_extension<NT,Root>& x , 
  const Sqrt_extension<NT,Root>& y) {
    return Min<Sqrt_extension<NT,Root> >() (x,y);
}

template< class NT, class Root >
inline Sqrt_extension<NT,Root> max BOOST_PREVENT_MACRO_SUBSTITUTION
( const Sqrt_extension<NT,Root>& x , 
  const Sqrt_extension<NT,Root>& y) {
    return Min<Sqrt_extension<NT,Root> >() (x,y);
}

namespace CGALi {


// Representation of the Bitstream Descartes traits class
template<typename Coefficient_,typename AlgebraicReal, 
         typename Integer_ >
class Bitstream_descartes_traits_on_vert_line_rep {
   
public:

    typedef AlgebraicReal Algebraic_real;
    typedef Integer_ Integer;
    typedef Coefficient_ Coefficient;

    // Coefficient is a polynomial type
    typedef typename Coefficient::NT Coeff_NT;

    typedef typename Algebraic_real::Rational Boundary;

    typedef typename Algebraic_real::Rational Rational;

    typedef typename CGAL::Coercion_traits<Coeff_NT,Rational>::Type
    Coercion;

    typedef typename CGAL::Coercion_traits<Coeff_NT,Rational>::Cast
    Rational_to_coercion_cast;

    typedef typename 
    CGAL::Get_arithmetic_kernel<Rational>::Arithmetic_kernel::
    Bigfloat_interval BFI;

    typedef typename 
    CGAL::Bigfloat_interval_traits<BFI>::Boundary BF;

    typedef CGAL::CGALi::Best_approximation_cache<Coefficient,Integer> 
    Best_approximation_cache;

    typedef Bitstream_descartes_traits_on_vert_line_rep
    <Coefficient,Algebraic_real,Integer> Self;

    Bitstream_descartes_traits_on_vert_line_rep
    (Algebraic_real alpha,bool use_approx_memory) 
	: alpha_(alpha), use_approx_memory(use_approx_memory)
    {	
    } 


    /*
      Bitstream_descartes_traits_on_vert_line_rep
      (const Self& s) 
      : alpha_(s.alpha)
      {	
      } 
    */

    bool uses_approx() const {
	return use_approx_memory;
    }

    Best_approximation_cache approx_mem() const {
	return approximations_;
    }
     

    Algebraic_real alpha() const {
	return alpha_;
    }

    class Approximator {

    private:
	
	Algebraic_real alpha_;

	Best_approximation_cache approximations_;

	bool use_approx_mem;

        long curr_prec;

    public:
	Approximator(Algebraic_real alpha) : alpha_(alpha),use_approx_mem(true), curr_prec(4) {};

	Approximator(Algebraic_real alpha,
		     Best_approximation_cache approximations) 
            : alpha_(alpha), approximations_(approximations), use_approx_mem(true), curr_prec(4) {};

	Approximator() {};

        template<typename NT>
        class Coeff_to_bfi_functor {

        public:

            typedef NT argument_type;
            typedef BFI result_type;

            BFI operator() (NT x) const {
                
                return CGAL::convert_to_bfi(x);

            }

        };
        

	typedef boost::numeric::interval<Coercion> Interval;

        long log_width(Algebraic_real alpha) const {
            CGAL_assertion(! alpha.is_rational());
            double d = CGAL::to_double(alpha.high()- alpha.low());
            if(d != 0.0) {
                return CGAL::CGALi::ceil_log2_abs(d);
            }
            return 0;
        }
    
#if CGAL_ACK_USE_NO_BFI_APPROX_IN_BITSTREAM_TRAITS
#warning uses no bfi-approx!

Integer operator() (Coefficient f, long p) {

            //AcX_DSTREAM("Called approximator with " << f << " and " << p << std::endl);
            //AcX_DSTREAM("Alpha: " << (CGAL::to_double(alpha_.high()-alpha_.low())) << " " << rational_to_double(alpha_.high()) << std::endl);


      
            if(use_approx_mem) {
                if(approximations_.approximation_exists(f)) {
                    long prec;
                    Integer approx;
                    approximations_.get_best_approximation(f,prec,approx);
                    if(prec>=p) {
                        //AcX_DSTREAM("use older approx.." << std::flush);
                        Integer divisor = CGAL::ipower((Integer)2,prec-p);
                        Integer return_value = CGAL::div(approx,divisor);
                        //AcX_DSTREAM("done" << std::endl);
                        return return_value;
                    }
                }
            }
            
            //AcX_DSTREAM("Have to approximate..." << std::endl);
	 
            Coercion bound;
            
            Rational_to_coercion_cast rational_to_coercion_cast;

            if(p<0) {
                bound = rational_to_coercion_cast
                    (Rational(CGAL::ipower((Integer)2,-p)));
            }
            else {
                typename CGAL::Fraction_traits<Rational>::Compose compose;
                bound =
                    rational_to_coercion_cast
                    (compose(1,CGAL::ipower(Integer(2),p)));
            }
            Interval f_alpha_iv;



            if(alpha_.is_rational()) {
                Coercion f_alpha=f.evaluate(alpha_.low());
                f_alpha_iv=Interval(f_alpha,f_alpha);
                
            }
            else {
                typedef Approximate_arithmetic_controller
                    <Coefficient,Algebraic_real>
                    Approximate_controller;

                Approximate_controller approx_controller(f,alpha_);

                while(true) {

                    f_alpha_iv = approx_controller.interval_approximation();
	      
                    //AcX_DSTREAM("done" << std::endl);
                    if(CGAL::compare
                       (f_alpha_iv.upper()-f_alpha_iv.lower(),bound) ==
                       CGAL::SMALLER) {
                        break;
                    }
                    //AcX_DSTREAM("refine..." << std::flush);
	      
                    approx_controller.refine_value();
	      
                    //AcX_DSTREAM("done" << std::endl);
                }
            }
       

            Integer return_value;
            typename CGAL::Fraction_traits<Coercion>::Decompose decompose;
            typedef typename 
                CGAL::Fraction_traits<Coercion>::Numerator_type Numerator;
            typedef typename 
                CGAL::Fraction_traits<Coercion>::Denominator_type Denominator;
            if(sign(f_alpha_iv.lower())==CGAL::ZERO) {
                return_value = Integer(0);
            }
            else if(CGAL::sign(f_alpha_iv.lower()) != 
                    CGAL::sign(f_alpha_iv.upper())) {
                return_value = Integer(0);
            }
            else {
                Numerator num;
                Denominator denom;

                decompose(f_alpha_iv.upper(),num,denom);                


                if(p<0) {
                    denom *= CGAL::ipower(Integer(2),-p);
                }
                else {
                    num *= CGAL::ipower(Integer(2),p);
                }
                return_value = CGAL::div(num,denom);

            }

            if(use_approx_mem) {
                approximations_.update_approximation(f,p,return_value);
            }
            //AcX_DSTREAM(" returning " << " " <<  return_value << std::endl);

            return return_value;
	}
    

#else

        CGAL::Polynomial<BFI> _convert_polynomial_to_bfi(Coefficient f) const {
            std::vector<BFI> coeffs;
            for(int i = 0; i <= f.degree(); i++) {
                coeffs.push_back(CGAL::convert_to_bfi(f[i]));
            }
            return CGAL::Polynomial<BFI>(coeffs.begin(), coeffs.end());   
        }
        
        
        Integer operator() (Coefficient f, long p) {
            
            typename CGAL::CGALi::Float_traits<BF>::Get_exponent get_exp;
            typename CGAL::CGALi::Float_traits<BF>::Get_mantissa get_m;

            long old_prec = CGAL::get_precision(BFI());
            
            CGAL::Polynomial<BFI> f_bfi;
            BFI alpha_bfi, f_alpha_bfi;
            
            long prec = 16;
            
            Integer return_value;
            
            long wbit = 0;

            while(true) {
                CGAL::set_precision(BFI(),prec);
                
                f_bfi = _convert_polynomial_to_bfi(f);
                alpha_bfi = CGAL::convert_to_bfi(alpha_);
                
                f_alpha_bfi = f_bfi.evaluate(alpha_bfi);
                
                if(!CGAL::singleton(f_alpha_bfi)) {
                    long ceil = CGAL::CGALi::ceil_log2_abs(f_alpha_bfi);
                    long signi = CGAL::get_significant_bits(f_alpha_bfi);
                    wbit   = ceil - signi + p;
                     
                } 
                
                if(wbit<-5 || CGAL::singleton(f_alpha_bfi)) {
                    break;
                } else {
                    prec*=2;
                }
            }
            BF lower = CGAL::lower(f_alpha_bfi);
            
            long shift = - (p + get_exp(lower)); 
            Integer bfi_m(get_m(lower)); 
             if( shift > 0 ){
                while(shift>63) {
                    bfi_m = (bfi_m >> 63);
                    shift-=63;
                }
                bfi_m = (bfi_m >> shift);
            }else{
                // add 0 bits 
                bfi_m = (bfi_m << -shift);   
            }   
            return_value = bfi_m;
            CGAL::set_precision(BFI(),old_prec);
            
            return return_value;
        }
#endif

    };


    class Lower_bound_log2_abs {

    private:
	Algebraic_real alpha_;

    public:
	Lower_bound_log2_abs(Algebraic_real alpha) : alpha_(alpha) {}
	Lower_bound_log2_abs() {};

	typedef boost::numeric::interval<Coercion> Interval;

	long operator() (Coefficient f) {
            //std::cout << "Called lower_bound_log2_abs with " 
            //          << f << std::flush;
          
            CGAL_assertion(! alpha_.is_root_of(f));

            Interval f_alpha_iv;
            Coercion abs_lower,abs_upper;

            typedef CGAL::CGALi::Approximate_arithmetic_controller
                <Coefficient,Algebraic_real> Approximate_controller;
          
            Approximate_controller approx_controller(f,alpha_);

            while(true) {
                f_alpha_iv = approx_controller.interval_approximation();
                CGAL::Sign lower_sign = CGAL::sign(f_alpha_iv.lower());
                if(CGAL::sign(f_alpha_iv.upper())==lower_sign) {
                    if(lower_sign==CGAL::POSITIVE) {
                        abs_lower=f_alpha_iv.lower();
                        abs_upper=f_alpha_iv.upper();
                    }
                    else {
                        abs_lower=CGAL::abs(f_alpha_iv.upper());
                        abs_upper=CGAL::abs(f_alpha_iv.lower());
                    }
    
                    BFI bfi_low = CGAL::convert_to_bfi(abs_lower),
                        bfi_high = CGAL::convert_to_bfi(abs_upper);
                    long lower_bound = CGAL::CGALi::floor_log2_abs(bfi_low),
                        upper_bound = CGAL::CGALi::ceil_log2_abs(bfi_high);
                    
                    CGAL_assertion(upper_bound>=lower_bound);
                    if(upper_bound-lower_bound <=2) {
                        //std::cout << "returning " << lower_bound << std::endl;
                        return lower_bound;
                    }
                }
                approx_controller.refine_value();
            }
	}

    };

    class Upper_bound_log2_abs_approximator {

    private:
	Algebraic_real alpha_;

        std::vector<long> coefficients_for_alpha;

	int refinements_of_alpha;

	bool zero_test_enabled;

	int refinement_limit;

	// Stores id of polynomials which are known to vanish (or not to 
	// vanish) at alpha
	std::vector<long> zeroes,non_zeroes;

    public:
	Upper_bound_log2_abs_approximator(Algebraic_real alpha) 
            : alpha_(alpha),refinements_of_alpha(0),zero_test_enabled(false)
        {}

        Upper_bound_log2_abs_approximator() {};

        typedef boost::numeric::interval<Boundary> Boundary_interval;
	typedef boost::numeric::interval<Coercion> Interval;

	bool initial_upper_bound
        (Coefficient f, long& ub_log2_abs,bool& is_certainly_zero) {
            return improve_upper_bound(f,ub_log2_abs,is_certainly_zero);
	}

	bool improve_upper_bound
        (const Coefficient f, long& ub_log2_abs,bool& is_certainly_zero) {
            //AcX_DSTREAM("improve upper bound.." << f << std::flush);
            if(std::find(coefficients_for_alpha.begin(),
                         coefficients_for_alpha.end(),
                         f.id())!=coefficients_for_alpha.end()) {
                coefficients_for_alpha.clear();
                //AcX_DSTREAM("refine.." << std::flush);
                alpha_.bisect();
                ++refinements_of_alpha;
                if(refinements_of_alpha>=refinement_limit) {
                    zero_test_enabled=true;
                }
                //AcX_DSTREAM("done.." << std::flush);
            }
            coefficients_for_alpha.push_back(f.id());
            if(zero_test_enabled) {
                if(std::find(zeroes.begin(),
                             zeroes.end(),
                             f.id())!=zeroes.end()) {
                    is_certainly_zero=true;
                    return true;
                }
                else if(std::find(non_zeroes.begin(),
                                  non_zeroes.end(),
                                  f.id())!=non_zeroes.end()) {
                    is_certainly_zero=false;
                }
                else {
                    bool zero = CGAL::CGALi::is_root_of(alpha_,f);
                    if(zero) {
                        zeroes.push_back(f.id());
                        is_certainly_zero=true;
                        return true;
                    }
                    else {
                        non_zeroes.push_back(f.id());
                        is_certainly_zero=false;
                    }
                }
            }
            else {
                is_certainly_zero=false;
            }
            Boundary_interval alpha_iv(alpha_.low(),alpha_.high());
            Interval f_alpha_iv = evaluate_iv(f,alpha_iv);
            CGAL::Sign lower_sign = CGAL::sign(f_alpha_iv.lower());
            bool iv_contains_zero = lower_sign 
                != CGAL::sign(f_alpha_iv.upper());
            Coercion abs_upper
                =(CGAL::abs(f_alpha_iv.lower())<CGAL::abs(f_alpha_iv.upper())) 
                ? CGAL::abs(f_alpha_iv.upper()) 
                : CGAL::abs(f_alpha_iv.lower());
            Coercion abs_lower(0);
            if(! iv_contains_zero) {
                abs_lower = (CGAL::abs(f_alpha_iv.lower())
                             < CGAL::abs(f_alpha_iv.upper())) 
                    ? CGAL::abs(f_alpha_iv.lower()) 
                    : CGAL::abs(f_alpha_iv.upper());
            }
            //Numerator upper_num;
            //Denominator upper_denom;
            //decompose(abs_upper,upper_num,upper_denom);
            if(CGAL::sign(abs_upper)==CGAL::ZERO) {
                is_certainly_zero=true;
                return true;
            }
            ub_log2_abs = CGAL::CGALi::ceil_log2_abs
                (CGAL::convert_to_bfi(abs_upper));
            //ub_log2_abs = num_ceil_log2_abs(upper_num) - 
            //    denom_floor_log2_abs(upper_denom);

            if(! iv_contains_zero) {
                //Numerator lower_num;
                //Denominator lower_denom;
                //decompose(abs_lower,lower_num,lower_denom);
                //long lb_log2_abs = num_floor_log2_abs(lower_num) 
                //    - denom_ceil_log2_abs(lower_denom);
                long lb_log2_abs 
                    = CGAL::CGALi::ceil_log2_abs
                    (CGAL::convert_to_bfi(abs_lower));
                CGAL_assertion(ub_log2_abs >= lb_log2_abs);
                //AcX_DSTREAM("Upper: " << ub_log2_abs << " Lower: " << lb_log2_abs << std::endl); 
                //AcX_DSTREAM(((ub_log2_abs - lb_log2_abs) <= 2) << std::endl);
	    
                return ((ub_log2_abs - lb_log2_abs) <= 2);
            }
            else {
                //AcX_DSTREAM("Upper: " << ub_log2_abs << std::endl);
                return false;
            }
	}
    };

private:
      
    Algebraic_real alpha_;

    Best_approximation_cache approximations_;

    bool use_approx_memory;

};

/*! 
 * \brief Traits class for the Bitstream Descartes method over algebraic 
 * extensions.
 *
 * Approximates coefficients of the polynomial 
 * <tt>g=f|<sub>x=a</sub></tt> 
 * using interval arithmetic. The \c AlgebraicReal object \c a is fixed
 * for the traits class. Each coefficient is represented by an univariate
 * integer polynomial. For a more detailled specification of its member 
 * functions, see the documentation of 
 * \c CGAL::CGALi::Bitstream_descartes_rndl_tree
 *
 * Using this traits class, bivariate polynomials can be isolated with 
 * AcX::Bitstream_descartes_bfs, the inner variable \c x is interpreted as
 * the algebraic number \c a.
 */
template<typename Coefficient_,
         typename AlgebraicReal,
         typename Integer_
             = typename CGAL::Get_arithmetic_kernel<typename CGAL::Polynomial_traits_d<Coefficient_>::Innermost_coefficient>
             ::Arithmetic_kernel::Integer>
class Bitstream_descartes_traits_on_vert_line
    : public ::CGAL::Handle_with_policy
    <CGAL::CGALi::Bitstream_descartes_traits_on_vert_line_rep
      <Coefficient_,AlgebraicReal,Integer_> > {
      
public:
      
    //! Coefficient type
    typedef Coefficient_ Coefficient;

    //! Algebraic real type
    typedef AlgebraicReal Algebraic_real;

    //! Integer type
    typedef Integer_ Integer;
      
    typedef Bitstream_descartes_traits_on_vert_line
    <Coefficient,Algebraic_real,Integer> Self;
      
    typedef CGAL::CGALi::Bitstream_descartes_traits_on_vert_line_rep
    <Coefficient,AlgebraicReal,Integer> Rep;
      
    typedef  ::CGAL::Handle_with_policy<Rep> Handle;

    typedef typename Rep::Boundary Boundary;
      
    //! Functor type for the approximation of coefficients
    typedef typename Rep::Approximator Approximator;
      
    //! Functor type for lower bounds on coefficients
    typedef typename Rep::Lower_bound_log2_abs Lower_bound_log2_abs;
      
    //! Functor type for upper bounds on coefficients
    typedef typename Rep::Upper_bound_log2_abs_approximator
    Upper_bound_log2_abs_approximator;
      
    typedef typename Rep::Rational Rational;
      
    // We can just take these from NT_traits
    typedef typename CGAL::Real_embeddable_traits<Integer>::Sign
    Sign;
    typedef typename CGAL::CGALi::Real_embeddable_extension<Integer>
        ::Ceil_log2_abs Ceil_log2_abs_Integer;
    typedef typename CGAL::CGALi::Real_embeddable_extension<long>
        ::Ceil_log2_abs Ceil_log2_abs_long;

      
    /*!
     * \brief Constructor for the polynomial <tt>f|<sub>x=a</sub></tt>
     *
     * Initialises an instance such that all bivariate polynomials are
     * interpreted as univariate polynomials in the algebraic extension
     * <tt>Z[a]</tt>.
     * If the \c use_approx_mem flag is set, the best known approximation
     * of each coefficient is stored in a 
     * \c CGAL::CGALi::Intern::Best_approximation_cache instance.
     */
    Bitstream_descartes_traits_on_vert_line
    (AlgebraicReal alpha=AlgebraicReal(),bool use_approx_mem = true) 
	: Handle(alpha,use_approx_mem)
    {	
    } 
      
    /*
      Bitstream_descartes_traits_on_vert_line
      (const Self& s) 
      : Handle(s)
      {	
      } 
    */
      
    Algebraic_real alpha() const {
        return this->ptr()->alpha();
    }

    //! See the documentation of CGAL::CGALi::Bitstream_descartes_rndl_tree.      
    Approximator approximator_object() const {
	if(this->ptr()->uses_approx()) {
            return Approximator(this->ptr()->alpha(),this->ptr()->approx_mem());
	}
	else {
            return Approximator(this->ptr()->alpha());
	}
    }
      
    //! See the documentation of CGAL::CGALi::Bitstream_descartes_rndl_tree.            
    Upper_bound_log2_abs_approximator 
    upper_bound_log2_abs_approximator_object() const {
	return Upper_bound_log2_abs_approximator(this->ptr()->alpha());
    }
      
    //! See the documentation of CGAL::CGALi::Bitstream_descartes_rndl_tree.            
    Lower_bound_log2_abs lower_bound_log2_abs_object() const {
	return Lower_bound_log2_abs(this->ptr()->alpha());
    }
      
      
    //! See the documentation of CGAL::CGALi::Bitstream_descartes_rndl_tree. 
    class Boundary_creator {
	
    public:
	
	Boundary_creator() {}
	
	Boundary operator() (Integer x,long p) {
            Integer num=x, denom,two(2),q,r;
            if(p <0) {
                CGAL::div_mod(num,two,q,r);
                while(r==Integer(0) && p<0) {
                    num=q;
                    p++;
                    CGAL::div_mod(num,two,q,r);
                }
                denom = CGAL::ipower(Integer(2),-p);
            }
            else {
                num*=CGAL::ipower(Integer(2),p);
                denom=1;
            }
            Rational b(num,denom);
            CGAL::simplify(b);
            return b;
	} 
    };
      
};
    
} // namespace CGALi

CGAL_END_NAMESPACE

#endif

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

    typedef typename Algebraic_real::Rational Boundary;

    typedef typename Algebraic_real::Rational Rational;

    typedef typename 
    CGAL::Get_arithmetic_kernel<Rational>::Arithmetic_kernel::
        Bigfloat_interval BFI;

    typedef typename 
        CGAL::Bigfloat_interval_traits<BFI>::Boundary BF;

    typedef Bitstream_descartes_traits_on_vert_line_rep
        <Coefficient,Algebraic_real,Integer> Self;

    Bitstream_descartes_traits_on_vert_line_rep
    (Algebraic_real alpha) 
	: alpha_(alpha)
    {	
    } 

    Algebraic_real alpha() const {
	return alpha_;
    }

private:

    CGAL::Polynomial<BFI> _convert_polynomial_to_bfi(Coefficient f) const {
            std::vector<BFI> coeffs;
            for(int i = 0; i <= f.degree(); i++) {
                coeffs.push_back(CGAL::convert_to_bfi(f[i]));
            }
            return CGAL::Polynomial<BFI>(coeffs.begin(), coeffs.end());   
        }

    BFI convert_coefficient_to_bfi(Coefficient f) const {
        CGAL::Polynomial<BFI> f_bfi;
        BFI alpha_bfi, f_alpha_bfi;
        
        long p = CGAL::get_precision(BFI());

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
        CGAL::set_precision(BFI(),p);
        return f_alpha_bfi;
    }

    bool coefficient_vanishes(Coefficient f) const {
        return alpha_.is_root_of(f);
    }

public:

    class Approximator {

    private:
	
	const Bitstream_descartes_traits_on_vert_line_rep* _m_traits;

    public:
	Approximator
            (const Bitstream_descartes_traits_on_vert_line_rep* traits) 
            : _m_traits(traits) {};

	Approximator() {};

        Integer operator() (Coefficient f, long p) {
            
            typename CGAL::CGALi::Float_traits<BF>::Get_exponent get_exp;
            typename CGAL::CGALi::Float_traits<BF>::Get_mantissa get_m;

            long old_prec = CGAL::get_precision(BFI());

            CGAL::set_precision(BFI(),p);
            
            BFI f_alpha_bfi = _m_traits->convert_coefficient_to_bfi(f);

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
            CGAL::set_precision(BFI(),old_prec);
            
            return bfi_m;
        }
    };


    class Lower_bound_log2_abs {

    private:
	const Bitstream_descartes_traits_on_vert_line_rep* _m_traits;

    public:
	Lower_bound_log2_abs
        (const Bitstream_descartes_traits_on_vert_line_rep* traits) 
            : _m_traits(traits) {}
	
        Lower_bound_log2_abs() {};

	long operator() (Coefficient f) {
            //std::cout << "Called lower_bound_log2_abs with " 
            //          << f << std::flush;
          
            CGAL_assertion(! _m_traits->coefficient_vanishes(f));

            long old_prec = CGAL::get_precision(BFI());
            long prec = 4;

            BFI f_alpha_iv;

            long result;
            while(true) {
                CGAL::set_precision(BFI(),prec);
                f_alpha_iv = _m_traits->convert_coefficient_to_bfi(f);
                CGAL::Sign lower_sign = CGAL::sign(CGAL::lower(f_alpha_iv));
                if(CGAL::sign(CGAL::upper(f_alpha_iv))==lower_sign) {
                    BF abs_lower, abs_upper;
                    if(lower_sign==CGAL::POSITIVE) {
                        abs_lower=CGAL::lower(f_alpha_iv);
                        abs_upper=CGAL::upper(f_alpha_iv);
                    }
                    else {
                        abs_lower=CGAL::abs(CGAL::upper(f_alpha_iv));
                        abs_upper=CGAL::abs(CGAL::upper(f_alpha_iv));
                    }
                    long lower_bound = CGAL::CGALi::floor_log2_abs(abs_lower),
                        upper_bound = CGAL::CGALi::ceil_log2_abs(abs_upper);
                    CGAL_assertion(upper_bound>=lower_bound);
                    if(upper_bound-lower_bound <=2) {
                        result = lower_bound;
                        break;
                    }
                }
                prec*=2;

            }

            //std::cout << "returning " << result 
            //          << std::endl;
            CGAL::set_precision(BFI(),old_prec);

            return result;
        }

    };

    class Upper_bound_log2_abs_approximator {

    private:
	const Bitstream_descartes_traits_on_vert_line_rep* _m_traits;

	// Stores id of polynomials which are known to vanish (or not to 
	// vanish) at alpha
	std::vector<long> zeroes,non_zeroes;

        std::vector<long> coeffs_for_alpha;

        // Stores the last known approximation to ensure an improvement
        long prec;

    public:
	Upper_bound_log2_abs_approximator
        (const Bitstream_descartes_traits_on_vert_line_rep* traits) 
            : _m_traits(traits), prec(4)
        {}

        Upper_bound_log2_abs_approximator() : prec(4) {};

        typedef boost::numeric::interval<Boundary> Boundary_interval;

	bool initial_upper_bound
        (Coefficient f, long& ub_log2_abs,bool& is_certainly_zero) {
            return improve_upper_bound(f,ub_log2_abs,is_certainly_zero);
	}

	bool improve_upper_bound
        (const Coefficient f, long& ub_log2_abs,bool& is_certainly_zero) {
            //std::cout << "improve upper bound.." 
            // <<  _m_traits->alpha().id() << ", " << f << std::endl;

            long old_prec = CGAL::get_precision(BFI());

            if(std::find(zeroes.begin(),
                         zeroes.end(),
                         f.id())!=zeroes.end()) {
                //std::cout << "ZERO FROM CACHE" << std::endl;
                is_certainly_zero=true;
                return true;
            }
            else if(std::find(non_zeroes.begin(),
                              non_zeroes.end(),
                              f.id())!=non_zeroes.end()) {
                //std::cout << "NON-ZERO FROM CACHE" << std::endl;
                is_certainly_zero=false;
            }
            else {
                bool zero = _m_traits->coefficient_vanishes(f);
                if(zero) {
                    //std::cout << "THAT IS ZERO!" << std::endl;
                    zeroes.push_back(f.id());
                    is_certainly_zero=true;
                    return true;
                }
                else {
                    //std::cout << "THAT IS NOT ZERO!" << std::endl;
                    non_zeroes.push_back(f.id());
                    is_certainly_zero=false;
                }
            }
            if(std::find(coeffs_for_alpha.begin(),
                         coeffs_for_alpha.end(),f.id())!=
               coeffs_for_alpha.end()) {
                prec*=2;
                coeffs_for_alpha.clear();
            }
            coeffs_for_alpha.push_back(f.id());
            
            BFI f_alpha_iv = _m_traits->convert_coefficient_to_bfi(f);
            
            BF abs_upper = std::max(CGAL::abs(CGAL::lower(f_alpha_iv)),
                                    CGAL::abs(CGAL::upper(f_alpha_iv)));
            
            if(CGAL::sign(abs_upper)==CGAL::ZERO) {
                is_certainly_zero=true;
                CGAL::set_precision(BFI(),old_prec);
                return true;
            }

            ub_log2_abs = CGAL::CGALi::ceil_log2_abs(abs_upper);

            if(! CGAL::zero_in(f_alpha_iv) ) {
                
                BF abs_lower = std::min(CGAL::abs(CGAL::lower(f_alpha_iv)),
                                        CGAL::abs(CGAL::upper(f_alpha_iv)));
                long lb_log2_abs 
                    = CGAL::CGALi::floor_log2_abs
                    (CGAL::convert_to_bfi(abs_lower));
                CGAL_assertion(ub_log2_abs >= lb_log2_abs);
                CGAL::set_precision(BFI(),old_prec);
                return ((ub_log2_abs - lb_log2_abs) <= 2);
            }
            else {
                //std::cout << "Upper: " << ub_log2_abs << std::endl;
                CGAL::set_precision(BFI(),old_prec);
                return false;
            }
        }
    };

private:
      
    mutable Algebraic_real alpha_;

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
    (AlgebraicReal alpha=AlgebraicReal())
	: Handle(alpha)
    {	
    } 
      
    Algebraic_real alpha() const {
        return this->ptr()->alpha();
    }

    //! See the documentation of CGAL::CGALi::Bitstream_descartes_rndl_tree.      
    Approximator approximator_object() const {
        return Approximator(this->ptr());
    }
      
    //! See the documentation of CGAL::CGALi::Bitstream_descartes_rndl_tree.            
    Upper_bound_log2_abs_approximator 
    upper_bound_log2_abs_approximator_object() const {
	return Upper_bound_log2_abs_approximator(this->ptr());
    }
      
    //! See the documentation of CGAL::CGALi::Bitstream_descartes_rndl_tree.            
    Lower_bound_log2_abs lower_bound_log2_abs_object() const {
	return Lower_bound_log2_abs(this->ptr());
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

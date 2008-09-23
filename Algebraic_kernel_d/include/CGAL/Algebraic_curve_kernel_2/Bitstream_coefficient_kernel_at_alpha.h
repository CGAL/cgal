// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber  <mkerber@mpi-inf.mpg.de>
//
// ============================================================================
#ifndef CGAL_BITSTREAM_COEFFICIENT_KERNEL_AT_ALPHA_H
#define CGAL_BITSTREAM_COEFFICIENT_KERNEL_AT_ALPHA_H 1

CGAL_BEGIN_NAMESPACE

#include <CGAL/basic.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/convert_to_bfi.h>

namespace CGALi {

template < typename Polynomial_1_, typename AlgebraicReal_1 >
class Bitstream_coefficient_kernel_at_alpha;

template < typename Polynomial_1_, typename AlgebraicReal_1 >
class Bitstream_coefficient_kernel_at_alpha_rep {

public:

    typedef Polynomial_1_ Polynomial_1;

    typedef AlgebraicReal_1 Algebraic_real_1;

    Bitstream_coefficient_kernel_at_alpha_rep() {}

    Bitstream_coefficient_kernel_at_alpha_rep(Algebraic_real_1 alpha)
        : _m_alpha(alpha) {}

    friend class Bitstream_coefficient_kernel_at_alpha
        <Polynomial_1,Algebraic_real_1>;

private:
    Algebraic_real_1 _m_alpha;

    
};

template < typename Polynomial_1_, typename AlgebraicReal_1 >
class Bitstream_coefficient_kernel_at_alpha
    : public CGAL::Handle_with_policy
        < CGAL::CGALi::Bitstream_coefficient_kernel_at_alpha_rep
            <Polynomial_1_, AlgebraicReal_1 > 
        >
{

public:

    //! \name typedefs
    //! @{

    typedef AlgebraicReal_1 Algebraic_real_1;
    
    typedef Polynomial_1_ Polynomial_1;

    typedef Polynomial_1 Coefficient;

    typedef typename 
    CGAL::Get_arithmetic_kernel<typename Coefficient::NT>::Arithmetic_kernel
        ::Bigfloat_interval Bigfloat_interval;

    typedef typename 
    CGAL::Get_arithmetic_kernel<typename Coefficient::NT>::Arithmetic_kernel
        ::Integer Integer;

    typedef typename 
    CGAL::Get_arithmetic_kernel<typename Coefficient::NT>::Arithmetic_kernel
        ::Rational Boundary;

    typedef CGAL::Handle_with_policy
        < CGAL::CGALi::Bitstream_coefficient_kernel_at_alpha_rep
            <Polynomial_1,Algebraic_real_1 > 
        > Handle;

    typedef Bitstream_coefficient_kernel_at_alpha<Coefficient,
                                                  Algebraic_real_1> Self;

    //! @}

    //! \name Constructors
    // !@{

    Bitstream_coefficient_kernel_at_alpha() {}

    Bitstream_coefficient_kernel_at_alpha(Algebraic_real_1 alpha) 
        : Handle(alpha) {}

    //@}

    //! \name Functors
    //! @{

    struct Is_zero : public std::unary_function<Coefficient,bool> {
        
        Is_zero(Algebraic_real_1 alpha) : _m_alpha(alpha) {}

        bool operator() (Coefficient f) const {
            return _m_alpha.is_root_of(f);
        }

    private:

        Algebraic_real_1 _m_alpha;

    };

    Is_zero is_zero_object() const {
        return Is_zero(this->ptr()->_m_alpha);
    }

    struct Convert_to_bfi 
        : public std::unary_function<Coefficient,Bigfloat_interval> {
        
        Convert_to_bfi(Algebraic_real_1 alpha) : _m_alpha(alpha) {}

        Bigfloat_interval operator() (Coefficient f) const {
            typename CGAL::Polynomial_traits_d<Coefficient>
                ::template Rebind<Bigfloat_interval,1>::Other::Type f_bfi;
            Bigfloat_interval alpha_bfi, f_alpha_bfi;
            
            long p = CGAL::get_precision(Bigfloat_interval());
            
            long prec = 16;
            
            long wbit = 0;
            
            while(true) {
                CGAL::set_precision(Bigfloat_interval(),prec);
                
                f_bfi = this->_convert_polynomial_to_bfi(f);
                alpha_bfi = CGAL::convert_to_bfi(_m_alpha);
                
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
            CGAL::set_precision(Bigfloat_interval(),p);
            return f_alpha_bfi;
        }
        
    private:
        
        typename CGAL::Polynomial_traits_d<Coefficient>
        ::template Rebind<Bigfloat_interval,1>::Other::Type
        _convert_polynomial_to_bfi(Coefficient f) const {

            typename
                CGAL::Polynomial_traits_d<Coefficient>::Degree degree;
            typename
                CGAL::Polynomial_traits_d<Coefficient>::Get_coefficient coeff;
            std::vector<Bigfloat_interval> coeffs;
            for(int i = 0; i <= degree(f); i++) {
                coeffs.push_back(CGAL::convert_to_bfi(coeff(f,i)));
            }
            return typename CGAL::Polynomial_traits_d<Coefficient>
                ::template Rebind<Bigfloat_interval,1>::Other
                ::Construct_polynomial()(coeffs.begin(),coeffs.end());   
        }
        
        Algebraic_real_1 _m_alpha;
    };

    Convert_to_bfi convert_to_bfi_object() const {
        return Convert_to_bfi(this->ptr()->_m_alpha);
    }

    // @}

};

} // namespace CGALi

CGAL_END_NAMESPACE


#endif // CGAL_BITSTREAM_COEFFICIENT_KERNEL_AT_ALPHA_H

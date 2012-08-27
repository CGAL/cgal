// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// Author(s)     : Michael Kerber  <mkerber@mpi-inf.mpg.de>
//
// ============================================================================
#ifndef CGAL_BITSTREAM_COEFFICIENT_KERNEL_AT_ALPHA_H
#define CGAL_BITSTREAM_COEFFICIENT_KERNEL_AT_ALPHA_H 1

namespace CGAL {

#include <CGAL/basic.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/convert_to_bfi.h>

namespace internal {

template < typename AlgebraicKernel_1 >
class Bitstream_coefficient_kernel_at_alpha;

template < typename AlgebraicKernel_1 >
class Bitstream_coefficient_kernel_at_alpha_rep {

public:

    typedef AlgebraicKernel_1 Algebraic_kernel_d_1;

    typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;

    typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;

    Bitstream_coefficient_kernel_at_alpha_rep() {}

    Bitstream_coefficient_kernel_at_alpha_rep(Algebraic_kernel_d_1* kernel,
                                              Algebraic_real_1 alpha)
        : _m_kernel(kernel), _m_alpha(alpha) {}

    friend class Bitstream_coefficient_kernel_at_alpha
        <Algebraic_kernel_d_1>;

private:
    Algebraic_kernel_d_1* _m_kernel;
    Algebraic_real_1 _m_alpha;

    
};

template < typename AlgebraicKernel_1 >
class Bitstream_coefficient_kernel_at_alpha
    : public CGAL::Handle_with_policy
        < CGAL::internal::Bitstream_coefficient_kernel_at_alpha_rep
            <AlgebraicKernel_1 > 
        >
{

public:

    //! \name typedefs
    //! @{

    typedef AlgebraicKernel_1 Algebraic_kernel_d_1;

    typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;
    
    typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;

    typedef Polynomial_1 Coefficient;

    typedef typename 
    CGAL::Get_arithmetic_kernel<typename Coefficient::NT>::Arithmetic_kernel
    Arithmetic_kernel; 

    typedef typename Arithmetic_kernel::Integer Integer;

    typedef typename Algebraic_kernel_d_1::Bound Bound;
  
    typedef typename Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;

    typedef CGAL::internal::Bitstream_coefficient_kernel_at_alpha_rep
            <Algebraic_kernel_d_1>                                      Rep;
    typedef CGAL::Handle_with_policy<Rep>                               Base;
    typedef Bitstream_coefficient_kernel_at_alpha<Algebraic_kernel_d_1> Self;

    //! @}

public:
    //! \name Constructors
    // !@{

    Bitstream_coefficient_kernel_at_alpha() : Base(Rep()) {}

    Bitstream_coefficient_kernel_at_alpha(const Self& traits)
      : Base(static_cast<const Base&>(traits)) {}

    Bitstream_coefficient_kernel_at_alpha(Algebraic_kernel_d_1* kernel,
                                          Algebraic_real_1 alpha) 
      : Base(kernel,alpha) {}

    //@}

    //! \name Functors
    //! @{

    struct Is_zero : public std::unary_function<Coefficient,bool> {
        
        Is_zero(Algebraic_kernel_d_1* kernel,Algebraic_real_1 alpha) 
            : _m_kernel(kernel),_m_alpha(alpha) {}

        bool operator() (Coefficient f) const {
            return _m_kernel->is_zero_at_1_object() (f,_m_alpha);
        }

    private:
        Algebraic_kernel_d_1* _m_kernel;
        Algebraic_real_1 _m_alpha;

    };

    Is_zero is_zero_object() const {
        return Is_zero(this->ptr()->_m_kernel,this->ptr()->_m_alpha);
    }

    struct Convert_to_bfi 
        : public std::unary_function<Coefficient,Bigfloat_interval> {
        
        Convert_to_bfi(Algebraic_kernel_d_1* kernel,
		       Algebraic_real_1 alpha) 
	  : _m_kernel(kernel), _m_alpha(alpha) {}

        Bigfloat_interval operator() (Coefficient f) const {
            typename CGAL::Polynomial_traits_d<Coefficient>
                ::template Rebind<Bigfloat_interval,1>::Other::Type f_bfi;
            
	    typename Algebraic_kernel_d_1::Approximate_relative_1 approx_alpha
	      =_m_kernel->approximate_relative_1_object();

	    typedef typename Algebraic_kernel_d_1::Bound Bound;
	    
	    Bigfloat_interval alpha_bfi, f_alpha_bfi;
            
            long p = CGAL::get_precision(Bigfloat_interval());
            
            long prec = 16;
            
            long wbit = 0;
	    
	    while(true) {
                CGAL::set_precision(Bigfloat_interval(),prec);
                
                f_bfi = this->_convert_polynomial_to_bfi(f);
                
		std::pair<Bound,Bound> alpha_bounds
		  = approx_alpha(_m_alpha,prec);
		
		alpha_bfi = CGAL::hull
      		              (CGAL::convert_to_bfi(alpha_bounds.first),
			       CGAL::convert_to_bfi(alpha_bounds.second));
		
		f_alpha_bfi = f_bfi.evaluate(alpha_bfi);
                
                if(!CGAL::singleton(f_alpha_bfi)) {
                    long ceil = CGAL::internal::ceil_log2_abs(f_alpha_bfi);
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
                CGAL::Polynomial_traits_d<Coefficient>::Get_coefficient coeff;
            std::vector<Bigfloat_interval> coeffs(CGAL::degree(f)+1);
            const int d = CGAL::degree(f);
            for(int i = 0; i <= d; i++) {
                coeffs[i] = CGAL::convert_to_bfi(coeff(f,i));
            }
            return typename CGAL::Polynomial_traits_d<Coefficient>
                ::template Rebind<Bigfloat_interval,1>::Other
                ::Construct_polynomial()(coeffs.begin(),coeffs.end());   
        }

        Algebraic_kernel_d_1* _m_kernel;
            
        Algebraic_real_1 _m_alpha;
    };

    Convert_to_bfi convert_to_bfi_object() const {
      return Convert_to_bfi(this->ptr()->_m_kernel,this->ptr()->_m_alpha);
    }

    // @}

};

} // namespace internal

} //namespace CGAL


#endif // CGAL_BITSTREAM_COEFFICIENT_KERNEL_AT_ALPHA_H

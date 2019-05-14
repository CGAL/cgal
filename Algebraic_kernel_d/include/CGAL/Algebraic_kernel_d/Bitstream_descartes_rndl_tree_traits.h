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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================
#ifndef CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_RNDL_TREE_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_RNDL_TREE_TRAITS_H

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/Bigfloat_interval_traits.h>
#include <CGAL/Algebraic_kernel_d/Float_traits.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/convert_to_bfi.h>
#include <CGAL/tss.h>

#include <vector>

#if CGAL_USE_CORE
namespace CORE { class BigInt; } 
#endif

namespace CGAL {

namespace internal {

#if CGAL_USE_CORE
// bugfix for CORE by Michael Kerber 
// why is there a specialized function for CORE?
inline CORE::BigInt shift_integer_by(CORE::BigInt x, long shift){
  if( shift > 0 ){
    while(shift>63) {
      x = (x >> 63);
      shift-=63;
    }  
    x = (x >> shift);
  }else{
    // add 0 bits 
    x = (x << -shift);   
  }   
  return x; 
}
#endif

template <class Shiftable>
Shiftable shift_integer_by(Shiftable x, long shift){
  if( shift > 0 ){
    x >>= shift;
  }else{
    x <<= -shift;  // adds 0 bits 
  }   
  return x; 
}

// forward
template <typename BitstreamCoefficientKernel> 
class Bitstream_descartes_rndl_tree_traits;


template <class BitstreamCoefficientKernel> 
class Bitstream_descartes_rndl_tree_traits_rep {

public:

    typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

    Bitstream_descartes_rndl_tree_traits_rep
        (Bitstream_coefficient_kernel kernel)
	: _m_kernel(kernel)
    {	
    } 

    Bitstream_descartes_rndl_tree_traits_rep() {} 

private:
    
    Bitstream_coefficient_kernel _m_kernel;

    friend class Bitstream_descartes_rndl_tree_traits
        <Bitstream_coefficient_kernel>;

}; // end of class Bitstream_descartes_rndl_tree_traits_rep

// A version that relies on a Bitstream_coefficient_kernel model
template <typename BitstreamCoefficientKernel>
class Bitstream_descartes_rndl_tree_traits
    : CGAL::Handle_with_policy
    <CGAL::internal::Bitstream_descartes_rndl_tree_traits_rep
        <BitstreamCoefficientKernel> >
{

public:
    //! typedefs
    //! @{
    
    typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;

    typedef typename Bitstream_coefficient_kernel::Coefficient Coefficient;

    typedef typename Bitstream_coefficient_kernel::Bigfloat_interval BFI;
    typedef typename CGAL::Bigfloat_interval_traits<BFI>::Bound BF; 
    
    typedef typename 
        CGAL::Polynomial_type_generator<Coefficient,1>::Type POLY; 
    typedef  Bitstream_descartes_rndl_tree_traits
        < Bitstream_coefficient_kernel > Self;

    typedef CGAL::Handle_with_policy
        <CGAL::internal::Bitstream_descartes_rndl_tree_traits_rep
            <Bitstream_coefficient_kernel> >
        Base;

    typedef typename Bitstream_coefficient_kernel::Integer  Integer; 
    typedef typename Bitstream_coefficient_kernel::Bound Bound;

    //! @}

private:
    static const Self& get_default_instance(){
      Bitstream_coefficient_kernel kernel;
      CGAL_STATIC_THREAD_LOCAL_VARIABLE(Self, x,kernel);
      return x;
    }

public:

    //! \name Constructors
    //! @{

    Bitstream_descartes_rndl_tree_traits(const Bitstream_coefficient_kernel& kernel)
      : Base(kernel){} 
 
    Bitstream_descartes_rndl_tree_traits()
      : Base(static_cast<const Base&>(get_default_instance())){}

    // explicit copy-constructor, required by VC9
    Bitstream_descartes_rndl_tree_traits(const Self& traits)
      : Base(static_cast<const Base&>(traits)){}
  
    //! @}

    class Approximator {
        
    private:
	
	Bitstream_coefficient_kernel _m_kernel;

    public:
	Approximator
            (const Bitstream_coefficient_kernel& kernel) 
            : _m_kernel(kernel) {};

	Approximator() {};

        Integer operator() (Coefficient f, long p) {

            //std::cout << "Called approximator with f=" << f
            //          << " and p=" << p << std::endl;
            
            typename CGAL::internal::Float_traits<BF>::Get_exponent get_exp;
            typename CGAL::internal::Float_traits<BF>::Get_mantissa get_m;

            long old_prec = CGAL::get_precision(BFI());
            long prec = 4;

            BFI f_alpha_bfi;
            
            while(true) {
                
                CGAL::set_precision(BFI(),prec);
                
                f_alpha_bfi = _m_kernel.convert_to_bfi_object()(f);
                if(CGAL::singleton(f_alpha_bfi)) {
                    break;
                }
                if(CGAL::internal::ceil_log2_abs(CGAL::upper(f_alpha_bfi)-
                                              CGAL::lower(f_alpha_bfi)) <=-p) {
                    break;
                } else {
                    prec*=2;
                }
                
            }

            BF lower = CGAL::lower(f_alpha_bfi);
                
            long shift = - (p + get_exp(lower)); 
            Integer bfi_m(get_m(lower)); 
            bfi_m = shift_integer_by(bfi_m,shift);
           
//             if( shift > 0 ){               
//               while(shift>63) { // this is a bug fix HACK for CORE::BigInt
//                  bfi_m = (bfi_m >> 63);
//                  shift-=63;
//                }
//                bfi_m = (bfi_m >> shift);
//             }else{
//                 // add 0 bits 
//                 bfi_m = (bfi_m << -shift);   
//             }  
            CGAL::set_precision(BFI(),old_prec);
            
            //std::cout << "returns " << bfi_m << std::endl;

            return bfi_m;
        }
    };

    Approximator approximator_object() const {
        return Approximator(this->ptr()->_m_kernel);
    }


    class Lower_bound_log2_abs {

    private:
	Bitstream_coefficient_kernel _m_kernel;

    public:
	Lower_bound_log2_abs
        (const Bitstream_coefficient_kernel& kernel) 
            : _m_kernel(kernel) {}
	
        Lower_bound_log2_abs() {};

	long operator() (Coefficient f) {
            //std::cout << "Called lower_bound_log2_abs with " 
            //          << f << std::flush;
          
            CGAL_assertion(! _m_kernel.is_zero_object()(f));

            long old_prec = CGAL::get_precision(BFI());
            long prec = 4;

            BFI f_alpha_iv;

            long result;
            while(true) {
                CGAL::set_precision(BFI(),prec);
                f_alpha_iv = _m_kernel.convert_to_bfi_object()(f);
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
                    long lower_bound = CGAL::internal::floor_log2_abs(abs_lower),
                        upper_bound = CGAL::internal::ceil_log2_abs(abs_upper);
                    CGAL_assertion(upper_bound>=lower_bound);
                    if(upper_bound-lower_bound <=2) {
                        result = lower_bound;
                        break;
                    }
                }
                prec*=2;

            }

            //std::cout << "returning " << result << std::endl;
            CGAL::set_precision(BFI(),old_prec);

            return result;
        }

    };
    
    Lower_bound_log2_abs lower_bound_log2_abs_object() const {
	return Lower_bound_log2_abs(this->ptr()->_m_kernel);
    }


    class Upper_bound_log2_abs_approximator {

    private:
	Bitstream_coefficient_kernel _m_kernel;

	// Stores id of polynomials which are known to vanish (or not to 
	// vanish) at alpha
	std::vector<Coefficient> zeroes,non_zeroes;

        std::vector<Coefficient> coeffs_for_alpha;

        // Stores the last known approximation to ensure an improvement
        long prec;

    public:
	Upper_bound_log2_abs_approximator
        (const Bitstream_coefficient_kernel& kernel) 
            : _m_kernel(kernel), prec(4)
        {}

        Upper_bound_log2_abs_approximator() : prec(4) {};

	bool initial_upper_bound
        (Coefficient f, long& ub_log2_abs,bool& is_certainly_zero) {
            return improve_upper_bound(f,ub_log2_abs,is_certainly_zero);
	}

	bool improve_upper_bound
        (const Coefficient f, long& ub_log2_abs,bool& is_certainly_zero) {
            //std::cout << "improve upper bound.." 
            // << f << std::endl;

            long old_prec = CGAL::get_precision(BFI());

            if(std::find(zeroes.begin(),
                         zeroes.end(),
                         f)!=zeroes.end()) {
                //std::cout << "ZERO FROM CACHE" << std::endl;
                is_certainly_zero=true;
                return true;
            }
            else if(std::find(non_zeroes.begin(),
                              non_zeroes.end(),
                              f)!=non_zeroes.end()) {
                //std::cout << "NON-ZERO FROM CACHE" << std::endl;
                is_certainly_zero=false;
            }
            else {
                bool zero = _m_kernel.is_zero_object()(f);
                if(zero) {
                    //std::cout << "THAT IS ZERO!" << std::endl;
                    zeroes.push_back(f);
                    is_certainly_zero=true;
                    return true;
                }
                else {
                    //std::cout << "THAT IS NOT ZERO!" << std::endl;
                    non_zeroes.push_back(f);
                    is_certainly_zero=false;
                }
            }
            if(std::find(coeffs_for_alpha.begin(),
                         coeffs_for_alpha.end(),f)!=
               coeffs_for_alpha.end()) {
                prec*=2;
                coeffs_for_alpha.clear();
            }
            coeffs_for_alpha.push_back(f);
            
            BFI f_alpha_iv = _m_kernel.convert_to_bfi_object()(f);
            
            BF abs_upper = (std::max)(CGAL::abs(CGAL::lower(f_alpha_iv)),
                                      CGAL::abs(CGAL::upper(f_alpha_iv)));
            
            if(CGAL::sign(abs_upper)==CGAL::ZERO) {
                is_certainly_zero=true;
                CGAL::set_precision(BFI(),old_prec);
                return true;
            }

            ub_log2_abs = CGAL::internal::ceil_log2_abs(abs_upper);

            if(! CGAL::zero_in(f_alpha_iv) ) {
                
              BF abs_lower = (std::min)(CGAL::abs(CGAL::lower(f_alpha_iv)),
                                        CGAL::abs(CGAL::upper(f_alpha_iv)));
                long lb_log2_abs 
                    = CGAL::internal::floor_log2_abs
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
    
    Upper_bound_log2_abs_approximator 
    upper_bound_log2_abs_approximator_object() const {
	return Upper_bound_log2_abs_approximator(this->ptr()->_m_kernel);
    }


    // TODO: Look whether this is best possible
    class Bound_creator {
	
    public:
	
	Bound_creator() {}
	
	Bound operator() (Integer x,long p) {
            Integer num=x, denom,two(2),q,r;
            if(p < 0) {
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
            Bound b(num);
	    b /= Bound(denom);
            CGAL::simplify(b);
            return b;
	} 
    };
    
    typedef typename CGAL::Real_embeddable_traits<Integer>::Sgn Sign;
    typedef typename CGAL::internal::Real_embeddable_extension<Integer>
        ::Ceil_log2_abs Ceil_log2_abs_Integer;
    typedef typename CGAL::internal::Real_embeddable_extension<long>
        ::Ceil_log2_abs Ceil_log2_abs_long;

}; // end of class Bitstream_descartes_rndl_tree_traits
    

} // namespace internal

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_RNDL_TREE_TRAITS_H

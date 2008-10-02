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

#ifndef TRIGONOMETRIC_APPROXIMATOR_H
#define TRIGONOMETRIC_APPROXIMATOR_H 1

#include <CGAL/basic.h>

#include <CGAL/Polynomial.h>
#include <CGAL/convert_to_bfi.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/ipower.h>
#include <CGAL/Algebraic_kernel_d/Float_traits.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi{

// So far, this is independent of the actual value to approximate
// Precondition: 0<x<Pi/4
template<typename Arithmetic_traits>
long number_of_summands_for_sin(long precision) {
    long n=0;
    typename Arithmetic_traits::Integer n_fac(1);
    while(CGAL::CGALi::floor_log2_abs(n_fac)<precision+1) {
        n++;
        n_fac*=n;
    }
    return n;
}

template<typename Arithmetic_traits>
long number_of_summands_for_arcsin(typename Arithmetic_traits::Rational s,
                                   long precision) {
    typedef typename Arithmetic_traits::Integer Integer;
    typedef typename Arithmetic_traits::Rational Rational;
    typedef typename Arithmetic_traits::Bigfloat_interval Bigfloat_interval;
    Rational s_pow = s*s*s;
    typename CGAL::Fraction_traits<Rational>::Compose compose;
    Rational bound 
        = compose(Integer(1),CGAL::ipower(Integer(2),precision+2))*
        (compose(1,1)-s*s);
    long m = 0;
    while(s_pow*compose(Integer(1),Integer(2*m+1))>bound) {
        m++;
        s_pow *= s*s;
    }
    return m;
}


} // namespace CGALi

template<typename Arithmetic_traits>
typename Arithmetic_traits::Bigfloat_interval pi(long precision)  {

    typedef typename Arithmetic_traits::Integer Integer;
    typedef typename Arithmetic_traits::Rational Rational;
    typedef typename Arithmetic_traits::Bigfloat_interval Bigfloat_interval;

    long old_prec = CGAL::get_precision(Bigfloat_interval());

    long error_offset = 2;
    Bigfloat_interval result;

    while(true) {
        
        long prec = 16;
        CGAL::set_precision(Bigfloat_interval(),prec);
        long m = 0;
        
        while((long)(precision-4*m-2) > 
              long(std::ceil(log(double(2*m+1))/log(2)))) {
            m++;
        }
        typedef typename CGAL::Bigfloat_interval_traits<Bigfloat_interval>
            ::Boundary Bigfloat_boundary;
        Bigfloat_boundary bb(Integer(1));
        bb = typename CGAL::CGALi::Float_traits<Bigfloat_boundary>
            ::Mul_by_pow_of_2() (bb,Integer(-precision-error_offset));
        Bigfloat_interval err = CGAL::hull(Bigfloat_interval(-bb), 
                                           Bigfloat_interval(bb));
        CGAL_assertion(CGAL::CGALi::ceil_log2_abs
                       (CGAL::upper(err)-CGAL::lower(err))==
                       -precision-error_offset+1);
                
        while(true) {
            Bigfloat_interval pow5 = CGAL::convert_to_bfi(Integer(1)),
                pow239 = CGAL::convert_to_bfi(Integer(1));
            Bigfloat_interval result1,result2;
            result = result1 = result2 = CGAL::convert_to_bfi(Integer(0));
            for(long i = 0; i < m; i++) {
                if(i%2==0) {
                    result1 += pow5/CGAL::convert_to_bfi(Integer(2*i+1));
                    result2 += pow239/CGAL::convert_to_bfi(Integer(2*i+1));
                } else {
                    result1 -= pow5/CGAL::convert_to_bfi(Integer(2*i+1));
                    result2 -= pow239/CGAL::convert_to_bfi(Integer(2*i+1));
                }
                pow5 /=CGAL::convert_to_bfi(Integer(25));
                pow239/=CGAL::convert_to_bfi(Integer(57121));
            }
            result=(result1*CGAL::convert_to_bfi(Integer(16)))
                / CGAL::convert_to_bfi(Integer(5))
                - (result2*CGAL::convert_to_bfi(Integer(4)))
                / CGAL::convert_to_bfi(Integer(239));
            // We don't expect that PI is rational...
            CGAL_assertion(! CGAL::singleton(result));
            
            if(CGAL::CGALi::ceil_log2_abs(CGAL::upper(result)-
                                          CGAL::lower(result))
               <=-precision-error_offset+1) { // precision -1 would suffices
                break;
            } else {
                prec*=2;
                CGAL::set_precision(Bigfloat_interval(),prec);
            }
        }
        result+=err;
        if(CGAL::CGALi::ceil_log2_abs(CGAL::upper(result)-
                                          CGAL::lower(result))
           <=-precision) {
            break;
        } else {
            error_offset++;
        }
    }
    CGAL::set_precision(Bigfloat_interval(),old_prec);
    return result;
    
}

template<typename Arithmetic_traits>
typename Arithmetic_traits::Bigfloat_interval arcsin
    (typename Arithmetic_traits::Rational x,long precision) {
    
    typedef typename Arithmetic_traits::Integer Integer;
    typedef typename Arithmetic_traits::Rational Rational;
    typedef typename Arithmetic_traits::Bigfloat_interval Bigfloat_interval;

    long old_prec = CGAL::get_precision(Bigfloat_interval());
    
    long error_offset = 2;
    Bigfloat_interval result;

    while(true) {
        
        long prec = 16;
        CGAL::set_precision(Bigfloat_interval(),prec);
        long m = CGAL::CGALi::number_of_summands_for_arcsin<Arithmetic_traits>
            (x,precision+error_offset);

        typedef typename CGAL::Bigfloat_interval_traits<Bigfloat_interval>
            ::Boundary Bigfloat_boundary;
        Bigfloat_boundary bb(Integer(1));
        bb = typename CGAL::CGALi::Float_traits<Bigfloat_boundary>
            ::Mul_by_pow_of_2() (bb,Integer(-precision-error_offset));
        Bigfloat_interval err = CGAL::hull(Bigfloat_interval(-bb), 
                                           Bigfloat_interval(bb));
        CGAL_assertion(CGAL::CGALi::ceil_log2_abs
                       (CGAL::upper(err)-CGAL::lower(err))==
                       -precision-error_offset+1);

        while(true) {
            CGAL::set_precision(Bigfloat_interval(),prec);
            Bigfloat_interval coeff,x_bfi,x_pow;
            coeff = CGAL::convert_to_bfi(Integer(1));
            x_bfi = CGAL::convert_to_bfi(x);
            x_pow = x_bfi;
            result = CGAL::convert_to_bfi(Integer(0));
            
            for(long n = 1; n < m; n++) {
                typename CGAL::Fraction_traits<Rational>::Compose compose;
                result+=x_pow*coeff;
                x_pow*=x_bfi*x_bfi;
                coeff=coeff*CGAL::convert_to_bfi(
                        compose(Integer((2*n-1))*Integer((2*n-1)),
                                Integer((2*n))*Integer((2*n+1))));
            }
            if(CGAL::singleton(result)) {
                break;
            }
            if(CGAL::CGALi::ceil_log2_abs(CGAL::upper(result)-
                                          CGAL::lower(result)) 
               < -(precision+1)) {
                break;
            } else {
                prec*=2;
            }
        }
        
        result+=err;
        if(CGAL::CGALi::ceil_log2_abs(CGAL::upper(result)-
                                      CGAL::lower(result)) 
           <= -precision) {
            break;
        } else {
            error_offset++;
        }
    }
    
    CGAL::set_precision(Bigfloat_interval(),old_prec);
    return result;
    
}

template<typename Arithmetic_traits, typename NT>
typename Arithmetic_traits::Bigfloat_interval sin
    (NT x,long precision) {
    
    typedef typename Arithmetic_traits::Integer Integer;
    typedef typename Arithmetic_traits::Rational Rational;
    typedef typename Arithmetic_traits::Bigfloat_interval Bigfloat_interval;

    long old_prec = CGAL::get_precision(Bigfloat_interval());
    
    long error_offset = 2;
    Bigfloat_interval result;

    while(true) {
        
        long prec = 16;
        CGAL::set_precision(Bigfloat_interval(),prec);
        long m = CGAL::CGALi::number_of_summands_for_sin<Arithmetic_traits>
            (precision+error_offset);

        typedef typename CGAL::Bigfloat_interval_traits<Bigfloat_interval>
            ::Boundary Bigfloat_boundary;
        Bigfloat_boundary bb(Integer(1));
        bb = typename CGAL::CGALi::Float_traits<Bigfloat_boundary>
            ::Mul_by_pow_of_2() (bb,Integer(-precision-error_offset));
        Bigfloat_interval err = CGAL::hull(Bigfloat_interval(-bb), 
                                           Bigfloat_interval(bb));
        CGAL_assertion(CGAL::CGALi::ceil_log2_abs
                       (CGAL::upper(err)-CGAL::lower(err))==
                       -precision-error_offset+1);

        
        while(true) {
            CGAL::set_precision(Bigfloat_interval(),prec);
            Bigfloat_interval coeff,x_bfi,x_pow;
            coeff = CGAL::convert_to_bfi(Integer(1));
            x_bfi = CGAL::convert_to_bfi(x);
            x_pow = x_bfi;
            result = CGAL::convert_to_bfi(Integer(0));
            
            for(long n = 1; n < m; n++) {
                typename CGAL::Fraction_traits<Rational>::Compose compose;
                result+=x_pow*coeff;
                x_pow*=x_bfi*x_bfi;
                coeff=coeff*CGAL::convert_to_bfi(compose(-1,(2*n)*(2*n+1)));
            }
            if(CGAL::singleton(result)) {
                break;
            }
            if(CGAL::CGALi::ceil_log2_abs(CGAL::upper(result)-
                                          CGAL::lower(result)) 
               < -precision-error_offset+1) {
                break;
            } else {
                prec*=2;
            }
        }
        result+=err;
        if(CGAL::CGALi::ceil_log2_abs(CGAL::upper(result)-
                                          CGAL::lower(result)) 
               <= -precision) {
            break;
        } else {
            error_offset++;
        }
    }
    CGAL::set_precision(Bigfloat_interval(),old_prec);
    return result;
    
}

template<typename Boundary>
std::pair<Boundary,Boundary>
approximate_sin_and_cos_of_angle(Boundary angle,long final_prec) {
    
    typedef typename 
            CGAL::Get_arithmetic_kernel<Boundary>::Arithmetic_kernel AT;
    typedef typename AT::Integer Integer;
    
    while(abs(angle)>180) {
        angle>0 ?  angle-=360 : angle+=360;
    }
    std::cout << angle << std::endl;
    bool between_90_and_270 = abs(angle)>=90;
    bool greater_180 =(angle<0);
    // Normalize
    angle=abs(angle);
    if(between_90_and_270) {
        angle=180-angle;
    }
    CGAL_assertion(angle>=0 && angle <= 90);
    
    bool greater_45 = (angle>=45);
    if(greater_45) {
        angle=90-angle;
    }
    
    Boundary sine, cosine;
    
    // Filter boundary case of 0 degree
    if(angle==Boundary(0) || 
       Boundary(1)/angle>CGAL::ipower(Integer(2),final_prec)) {
        sine = 0;
        cosine = 1;
    } else {
        
        typedef typename AT::Bigfloat_interval Bigfloat_interval;
        typedef typename Bigfloat_interval_traits<Bigfloat_interval>
            ::Boundary Bigfloat_boundary;
        
        
        long old_prec = CGAL::get_precision(Bigfloat_interval());
        
        long prec = 16;
        Boundary t;
        while(true) {
            CGAL::set_precision(Bigfloat_interval(),prec);
            std::cout << "increased prec to " << (prec) 
                      << std::endl; 
            Bigfloat_interval pi = CGAL::pi<AT>(prec);
            Bigfloat_interval s 
                = CGAL::sin<AT>
                (CGAL::median(
                         pi*CGAL::convert_to_bfi(angle)/
                         CGAL::convert_to_bfi(Integer(180))),
                 prec);
            Bigfloat_boundary x 
                = CGAL::median(CGAL::convert_to_bfi(Integer(1))/s + 
                               CGAL::sqrt
                               (CGAL::convert_to_bfi(Integer(1))/(s*s)-CGAL::convert_to_bfi(Integer(1))));
            
            int n = 0;
            typename CGAL::Fraction_traits<Boundary>::Compose compose;
            
            bool success=false;
            
            Bigfloat_boundary e0=x, e1=-1, olde0;
            Integer p0=0, q0=1, p1=1, q1=0,oldp0,oldq0;
            while(true) {
                Integer r = CGAL::CGALi::floor(e0/e1);
                Integer oldp0=p0, oldq0=q0;
                Bigfloat_boundary olde0=e0;
                e0=e1;p0=p1;q0=q1;
                e1=olde0-r*e1;
                p1=oldp0-r*p1;
                q1=oldq0-r*q1;
                if(q1!=Integer(0)) {
                    t = compose(p1,q1);
                    CGAL::simplify(t);
                    sine = CGAL::abs(2/(t+1/t));
                    CGAL::simplify(sine);
                    Bigfloat_interval asin 
                        = CGAL::arcsin<AT>(sine,prec)
                        * CGAL::convert_to_bfi(Integer(180))/pi;
                    
                    long bound = CGAL::CGALi::ceil_log2_abs
                        (CGAL::abs(asin-CGAL::convert_to_bfi(angle)));
                    success = (bound <= -final_prec);
                    typename 
                        CGAL::Coercion_traits<Boundary, 
                        Bigfloat_boundary>::Cast
                        cast;
                    if((cast(t)==cast(x)) || success) {
                        break;
                    }
                }
                n++;
            }
            if(success) {
                CGAL::set_precision(Bigfloat_interval(),old_prec);
                break;
            } else {
                prec*=2;
            }
        }
        
        cosine = (t-1/t)/(t+1/t);
        CGAL::simplify(cosine);
        
    }
    if(greater_45) {
        std::swap(sine,cosine);
    }
    
    if(greater_180) {
        sine = -sine;
    }
    if(between_90_and_270) {
        cosine=-cosine;
    }
    
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "sine=" << sine << std::endl;
    CGAL_ACK_DEBUG_PRINT << "cosine=" << cosine << std::endl;
#endif

    return std::make_pair(sine,cosine);


}

CGAL_END_NAMESPACE


#endif 

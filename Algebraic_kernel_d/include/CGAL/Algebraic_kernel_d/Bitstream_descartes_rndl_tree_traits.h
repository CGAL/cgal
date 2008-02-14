// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file NiX/Bitstream_descartes_rndl_tree_traits.h
  \brief Definition of \c NiX::Bitstream_descartes_rndl_tree_traits.
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_RNDL_TREE_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_RNDL_TREE_TRAITS_H

#define CGAL_INTERN_USE_BFI

#include <CGAL/basic.h>
#include <CGAL/Sqrt_extension.h>
// #include <CGAL/Algebraic_kernel_d/interval_support.h>
#include <CGAL/interval_support.h>

#include <CGAL/Algebraic_kernel_d/Integer_iterator.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/Algebraic_kernel_d/Float_traits.h>

#include <vector>

/*#ifdef LiS_HAVE_CORE
#include <NiX/CORE/BigInt.h>
#include <NiX/CORE/BigRat.h>
#include <NiX/CORE/Expr.h>
#endif // LiS_HAVE_CORE
*/

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template < class NT > class Bitstream_descartes_rndl_tree_traits;

// TODO: write a generic version that 
template <class NT>
class Bitstream_descartes_rndl_tree_traits{

// typedefs
    typedef typename CGALi::Get_arithmetic_kernel<NT>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat          BF;
    typedef typename AT::Bigfloat_interval BFI;
    typedef typename AT::Field_with_sqrt   FWS;
    
    typedef  CGAL::Polynomial<NT> POLY; 
    typedef  Bitstream_descartes_rndl_tree_traits< NT > Self;


template <class Coeff, class OI> 
inline 
void convert_coeffs(const CGAL::Polynomial<Coeff>& poly ,OI it) const {
    for(int i = 0; i <= poly.degree(); i++){
        *it++ = convert_to_bfi(poly[i]);
    }
}

template <class COEFF, class ROOT, class OI > 
inline  
void
convert_coeffs(
        const CGAL::Polynomial< Sqrt_extension<COEFF,ROOT> >& poly,
        OI it ) const {
    
    BFI root(0);
    for(int i = 0; i <= poly.degree(); i++){
        if(poly[i].is_extended()){
            //typename Coercion_traits<FWS,ROOT >::Cast cast_root;  
            //root = convert_to_bfi(NiX::sqrt(cast_root(poly[i].root())));
            root = CGALi::sqrt(convert_to_bfi(poly[i].root()));
            break;  
        }
    }
    
    for(int i = 0; i <= poly.degree(); i++){
        if(poly[i].is_extended()){
            *it++ = 
                convert_to_bfi(poly[i].a0()) + 
                convert_to_bfi(poly[i].a1()) * root ;
        }else{
            *it++ = convert_to_bfi(poly[i].a0());
        }
    } 
}

public:
    typedef int                   Coefficient;  
    typedef typename AT::Integer  Integer; 
    typedef typename AT::Rational Boundary;

private:   
//member
    POLY polynomial;
    mutable std::vector<BFI>          polynomial_approx; 
    mutable long                      current_prec; 
    
    void refine_approximation() const {
        CGAL_precondition(current_prec > 1);
        current_prec *= 2;
        // std::cout <<"ALGREAL: refine approx: "<<  current_prec<<std::endl;
        long old_prec = set_precision(BF(),current_prec);
        polynomial_approx.clear();
        convert_coeffs(
                polynomial,
                std::back_inserter(polynomial_approx));
        set_precision(BF(),old_prec);
    };
    
public:
    typedef CGALi::Integer_iterator<int> IIterator; 
    IIterator begin(){return IIterator(0);};
    IIterator end()  {return IIterator(polynomial_approx.size());}

    Bitstream_descartes_rndl_tree_traits(const POLY& p):polynomial(p){
        current_prec = 60; 
        long old_prec = set_precision(BF(), current_prec); 
        convert_coeffs(polynomial,std::back_inserter(polynomial_approx));
        set_precision(BF(),old_prec);
    };
    
    struct Boundary_creator
    //: public CGAL::Binary_function<Integer,long,Boundary>
    {
        typedef Integer first_argument_type;
        typedef long second_argument_type;
        typedef Boundary result_type; 
        Boundary operator()(Integer i, long e){ 
            Boundary power = Integer(1) << CGAL::abs(e);
            if( e < 0 ) power = 1 / power; 
            return Boundary(i) * power; 
        }
    };
      
    struct Lower_bound_log2_abs{             

      const Self* ptr;
      
      typedef int argument_type;
      typedef long result_type;
      Lower_bound_log2_abs(const Self* ptr_):ptr(ptr_){};
      result_type operator() (int i) {             
        CGAL_precondition(ptr->polynomial[i] != NT(0));
        while(CGALi::in_zero(ptr->polynomial_approx[i])){
          ptr->refine_approximation();
        }             
        typename CGALi::Real_embeddable_extension<BFI>::Floor_log2_abs floor_log2_abs;
        return floor_log2_abs(ptr->polynomial_approx[i]);
      }
    }; 
    
    struct Upper_bound_log2_abs{
        const Self* ptr; 
        
        Upper_bound_log2_abs(const Self* ptr_):ptr(ptr_){};
        
        bool initial_upper_bound(int i, long& upper_log, bool& is_certainly_zero){
            if (is_certainly_zero = ( ptr->polynomial[i] == NT(0) )) return true;
 
            typename CGALi::Real_embeddable_extension<BFI>::Ceil_log2_abs ceil_log2_abs;
            upper_log = ceil_log2_abs(ptr->polynomial_approx[i]);
            return true;
        }
        
        bool improve_upper_bound(int, long&, bool&){
            CGAL_precondition_msg(false,"this call is not needed");
            return true;
        }
    };



    class Approximator {
    public:
        Approximator(const Self* ptr_):ptr(ptr_){};
        const Self * ptr; 

        Integer operator() (int i, long p) { 
#if 1
          //std::cout << " Approximator begin  " << "Precision: " << p << std::flush;
            typename CGALi::Real_embeddable_extension<BFI>::Ceil_log2_abs ceil_log2_abs;
            typename CGALi::Float_traits<BF>::Get_exponent get_exp;
            typename CGALi::Float_traits<BF>::Get_mantissa get_m;
            
            BFI approx  = ptr->polynomial_approx[i];
            // if coeff = 0 -> fast return 
            if (CGALi::singleton(approx) && approx == BFI(0) ) {
              //std::cout << " Approximator end  " << std::endl;
              return Integer(0);
            }
            //get position of first wrong bit
            long wbit   = ceil_log2_abs(approx) - CGALi::get_significant_bits(approx)+p;
             // approx until pos of first wrong bit is negative
            while( wbit >= -5 && ! CGALi::singleton(approx) ){

                ptr->refine_approximation();
                approx = ptr->polynomial_approx[i]; 
                wbit   = ceil_log2_abs(approx) - CGALi::get_significant_bits(approx) + p;
            }  
            BF lower = CGALi::lower(approx); // could take upper also 
            long shift = - (p + get_exp(lower)); 
            Integer m(get_m(lower)); 
            if( shift > 0 ){
              while(shift>63) {
                m = (m >> 63);
                shift-=63;
              }
              m = (m >> shift);
            }else{
                // add 0 bits 
                CGAL_precondition(CGALi::singleton(approx));
                m = (m << -shift);   
            }     
            // std::cout << " Approximator end  " << std::endl;
 
            return m ;
#else
            //std::cout << " Approximator begin  " << std::endl;

            BFI approx = ptr->polynomial_approx[i];
            long shift = - (p + approx.exp()*14); 
            
            // std::cout <<"shift" <<  shift << std::endl;
            // shift must remove the error in m .. if there is one 
            while( shift  <= CORE::bitLength( approx.err())+3  
                    && !approx.isExact() ){
                //std::cout << prec << "," << approx.err() << std::endl;
                ptr->refine_approximation();
                approx = ptr->polynomial_approx[i]; 
                shift = - (p + approx.exp()*14); 
            }  
            Integer m; 
            if( shift > 0 ){
                m = approx.m() >> shift;   
            }else{
                // add 0 bits 
                CGAL_precondition(approx.isExact());
                m = approx.m() << -shift;   
            }     
            //std::cout << " Approximator end  " << std::endl;
 
            return m ;
#endif

        }
    };
    
   
    Lower_bound_log2_abs 
    lower_bound_log2_abs_object() const {
        return Lower_bound_log2_abs(this);
    }
    
    Upper_bound_log2_abs 
    upper_bound_log2_abs_object() const {
        return Upper_bound_log2_abs(this);
    }
    
    Approximator
    approximator_object() const {
        return Approximator(this);
    }

    typedef typename CGAL::Real_embeddable_traits<Integer>::Sign Sign;

    typedef typename CGALi::Real_embeddable_extension<Integer>::Ceil_log2_abs Ceil_log2_abs_Integer;
    typedef typename CGALi::Real_embeddable_extension<long>::Ceil_log2_abs Ceil_log2_abs_long;
};






} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_RNDL_TREE_TRAITS_H

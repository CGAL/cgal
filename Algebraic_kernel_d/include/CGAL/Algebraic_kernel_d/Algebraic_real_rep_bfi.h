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

// This is code is expreimental ! 

/*! \file NiX/Algebraic_real_rep_bfi.h
  \brief Algebraic_real_rep with refinement via interval rep 
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_REP_BFI_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_REP_BFI_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep.h>
/*#include <NiX/NT_traits.h>
#include <NiX/univariate_polynomial_utils.h>
#include <NiX/Algebraic_real_rep.h>
#include <NiX/Arithmetic_traits.h>
#include <NiX/interval_support.h>*/


CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class NT> class Get_arithmetic_kernel;

// definition of the Algebraic_real_rep_bfi x:
    
//IS_GENERAL: 
// low_  lower bound of x 
// high_ upper bound of x 
// polynomial_ a square free polynomial 
// sign_at_low_ = polynomial_.evaluate(low_)
// x is the only root of polynomial_ in the open interval ]low_,high_[ 
// low_ != x != high
// ******************* EXEPTION *******************
// x is rational: in this case low=high=x

 
template< class Coefficient_, class Field_> 
class Algebraic_real_rep_bfi
    : public Algebraic_real_rep<Coefficient_,Field_> {
    
    typedef Algebraic_real_rep<Coefficient_,Field_>  Base; 

    typedef Coefficient_                            Coefficient;
    typedef Field_                                  Field;
    
    typedef typename CGALi::Get_arithmetic_kernel<Coefficient>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat          BF;
    typedef typename AT::Bigfloat_interval BFI;
    typedef typename AT::Field_with_sqrt   FWS; 
    
    typedef CGAL::Polynomial<Coefficient>            Poly;
    typedef typename Poly::const_iterator           PIterator; 
    
    typedef Algebraic_real_rep_bfi <Coefficient,Field>     Self;   

    mutable std::vector<BFI>          polynomial_approx; 
    mutable long                      current_prec;

private:
    template <class NT, class OI> 
    inline 
    void convert_coeffs(const CGAL::Polynomial<NT>& poly ,OI it) const {
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
//                typename Coercion_traits<FWS,ROOT >::Cast cast_root;  
//                root = convert_to_bfi(NiX::sqrt(cast_root(poly[i].root())));
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
    //! creates the algebraic real from int \a i.
    explicit Algebraic_real_rep_bfi(int i = 0)
        :Base(i){}
 
    //! creates the algebraic real from Field \a m.
    explicit Algebraic_real_rep_bfi(const Field& m)
        :Base(m), current_prec(53) {}
    
    /*! \brief creates the algebraic real as the unique root of \a P
        in the open interval <var>]low,high[</var>.
        \pre the polynomial \a P is square free 
        \pre P(low)!=0
        \pre P(high)<0
        \pre x is the one and only root in the open interval of \a P.
    */
    Algebraic_real_rep_bfi(const Poly& P, Field LOW, Field HIGH):
        Base(P,LOW,HIGH), current_prec(53){};


    //! copy constructor
    Algebraic_real_rep_bfi(const Self& y)
        : Base(y), current_prec(53){}
    
    // assignment 
    Algebraic_real_rep_bfi& operator=(const Self& y) {
        if ( this != & y) {
            this->erase_from_list();
            this->copy_all_members(y);
            this->introduce(y);
            current_prec = y->current_prec;
        }
        NiX_expensive_postcond(this->self_test());
        return *this;
    }
    
   
    BFI evaluate_polynomial_approx(const BFI& x) const {
       // std::cout << "eval approx  begin"<< std::endl;
        typedef std::vector<BFI> BFI_VEC;
        typedef typename BFI_VEC::reverse_iterator RIT; 
        
        BFI result(0);
        for(RIT rit = polynomial_approx.rbegin();
            rit != polynomial_approx.rend();
            rit++){
            result = result * x + (*rit); 
        }
       // std::cout << "eval approx  end"<< std::endl;
        return result; 
    }
  
    void refine_approximation() const {
        CGAL_precondition(current_prec > 1);
        current_prec *= 2;
        // std::cout <<"ALGREAL: refine approx: "<<  current_prec<<std::endl;
        set_precision(BF(),current_prec);
        polynomial_approx.clear();
        convert_coeffs(
                this->polynomial(),
                std::back_inserter(polynomial_approx));
        
        Self const *next = static_cast<const Self *>(this->next);
        while(this != next){
            // std::cout << this << " " << next << std::endl;
            next->polynomial_approx = polynomial_approx; 
            next->current_prec      = current_prec; 
            next = static_cast<const Self *>(next->next);
        }
    };
    
    void refine() const{              
        if(this->is_rational()) return;
        
       // std::cout << "refine begin ------- "<< std::endl;

        long old_prec = set_precision(BF(),current_prec);
        
        if ( (int) polynomial_approx.size() != this->polynomial_.degree()+1) {
            refine_approximation();
        }        
        
        CGAL_postcondition(polynomial_approx.size() > 0);

        Field m = (this->low()+this->high())/Field(2);
        CGAL::simplify(m);

        BFI eval = evaluate_polynomial_approx(convert_to_bfi(m));
        //std::cout <<"approx eval: "<< eval << std::endl;
    
        //std::cout <<"true eval: "<< convert_to_bfi(this->polynomial().evaluate(m)) << std::endl;
        CGAL::Sign s = CGAL::sign(CGALi::lower(eval));
       
        // correct sign if needed
        if( s*CGAL::sign(CGALi::upper(eval) ) != CGAL::POSITIVE ){
            //std::cout << "APPROX FAILED-------------------------------"<<std::endl;
            s = this->polynomial().sign_at(m);
            if ( s != CGAL::ZERO ) {
                refine_approximation(); 
            }
        }
        
        CGAL_postcondition(s == this->polynomial_.sign_at(m));
       
        if ( s == CGAL::ZERO ) {
            this->learn_from(m);
        }else{
            if ( s == this->sign_at_low() ) this->low_  = m;
            else                            this->high_ = m; 
        }
        
        set_precision(BF(),old_prec);
        //std::cout << "refine end ----------------- "<<current_prec<<  std::endl;
    }
};

}//namespace CGALi

CGAL_END_NAMESPACE

#endif //CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_REP_BFI_H

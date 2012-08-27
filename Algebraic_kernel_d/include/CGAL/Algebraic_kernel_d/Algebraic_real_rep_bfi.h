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
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>
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
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep.h>
#include <CGAL/Interval_traits.h>
#include <CGAL/Bigfloat_interval_traits.h>
#include <CGAL/convert_to_bfi.h>
#include <CGAL/Sqrt_extension_fwd.h>
#include <CGAL/Get_arithmetic_kernel.h>

namespace CGAL {

// it would be nice to remove the explicit use of Sqrt_extension 
// in this file. However, it is much more efficient to convert the root 
// only once. see convert_to_bfi  
//template <class COEFF, class ROOT> class Sqrt_extension;  //use forward declaration instead

namespace internal {

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
    
    typedef typename CGAL::Get_arithmetic_kernel<Coefficient>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat          BF;
    typedef typename AT::Bigfloat_interval BFI;
    typedef typename AT::Field_with_sqrt   FWS; 
    
    typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type Poly;
    typedef typename Poly::const_iterator           PIterator; 
    
    typedef Algebraic_real_rep_bfi <Coefficient,Field>     Self;   

    mutable std::vector<BFI>          polynomial_approx; 
    mutable long                      current_prec;
        
private:
    template <class Polynomial_1, class OI> 
    inline 
    void convert_coeffs(const Polynomial_1& poly, OI it) const {
        typename CGAL::Polynomial_traits_d<Polynomial_1>::Get_coefficient
            coeff;
        for(int i = 0; i <= CGAL::degree(poly); i++){
            *it++ = convert_to_bfi(coeff(poly,i));
        }
    }

   template <class COEFF, class ROOT, class ACDE_TAG, class FP_TAG, class OI > 
    inline 
    void
    convert_coeffs(
            const CGAL::Polynomial< CGAL::Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG> >& poly,
            OI it ) const {
        
        BFI root(0);
        for(int i = 0; i <= CGAL::degree(poly); i++){
            if(poly[i].is_extended()){
//                typename Coercion_traits<FWS,ROOT >::Cast cast_root;  
//                root = convert_to_bfi(NiX::sqrt(cast_root(poly[i].root())));
                root = CGAL::sqrt(convert_to_bfi(poly[i].root()));
                break;  
            }
        }
        
        for(int i = 0; i <= CGAL::degree(poly); i++){
            if(poly[i].is_extended()){
                *it++ = 
                    convert_to_bfi(poly[i].a0()) + 
                    convert_to_bfi(poly[i].a1()) * root ;
            }else{
                *it++ = convert_to_bfi(poly[i].a0());
            }
        } 
    }

    typedef CGAL::Sign                              TRI_BOOL;

public:
    const Field&           low()           const {return this->low_;}
    const Field&           high()          const {return this->high_;}
    const Poly&            polynomial()    const {return this->polynomial_;}
    const TRI_BOOL&        sign_at_low()   const {return this->sign_at_low_;}
    
    template <class NTX>
    void strong_refine(const NTX& m) const{
        if(is_rational()) return;

        if( NTX(low()) <= m && m <= NTX(high()) ){
            CGAL_precondition(polynomial().sign_at(m)!=CGAL::ZERO);
            refine();
            while( NTX(low()) <= m && m <= NTX(high())) refine();
        }
    }

    bool                   is_rational()   const {return this->is_rational_;}
    const Field&           rational()      const {
        CGAL_precondition(is_rational());
        CGAL_precondition(low()==high());
        CGAL_precondition(polynomial().sign_at(low())==CGAL::ZERO);
        return this->low_;
    }

    CGAL::Comparison_result
    compare(const Field& y, bool are_distinct = false) const {  
        if(is_rational()) return CGAL::compare(rational(),y);
        if(y <= low()) {
            strong_refine(y);
            return CGAL::LARGER;
        }
        if(high() <= y) { 
            strong_refine(y);
            return CGAL::SMALLER;
        }
        
        // now: low < y < high 
        if(!are_distinct){
            if(sign_of_polynomial_at(y)==CGAL::ZERO){
                this->learn_from(y);
                return CGAL::EQUAL ;
            }
        }else{
            CGAL_precondition(polynomial().sign_at(y)!=CGAL::ZERO);
        }
        strong_refine(y);
        CGAL_postcondition(y < low() || high() < y );
        if(y < low()) return CGAL::LARGER;
        else          return CGAL::SMALLER;
    }
    
    CGAL::Comparison_result
    compare (const Algebraic_real_rep_bfi& y, bool are_distinct = false) const{
        if(  is_rational()) return -y.compare(  rational());
        if(y.is_rational()) return    compare(y.rational());
        
        // type of both x and y IS_GENERAL
        if ( high() <= y.low() ) return CGAL::SMALLER;
        if ( low()  >= y.high()) return CGAL::LARGER;

        // intersection isolating intervals is ]L,R[
        Field L = (low()  > y.low() ) ? low()  : y.low()  ;
        Field R = (high() < y.high()) ? high() : y.high() ;

        // refine to smaller intervals at intersection interval boundaries
        // this can change type() only to IS_RATIONAL
        this->refine_at(L);
        this->refine_at(R);
        y.refine_at(L);
        y.refine_at(R);

        if (  is_rational()) return -y.compare(  rational());
        if (y.is_rational()) return    compare(y.rational());

        // type of both x and y still IS_GENERAL
        if ( high()   <= y.low() ) return CGAL::SMALLER;
        if ( y.high() <=   low() ) return CGAL::LARGER;
        
        // filter 1 (optional): determine distinctness by refining intervals
#if NiX_REFINEMENTS_BEFORE_GCD > 0
        if (!are_distinct) {
            // we may want to refine a bit and hope for the best
            // because computing the gcd is expensive
            for (int ntries=0; ntries < NiX_REFINEMENTS_BEFORE_GCD; ntries++) {
                y.refine();
                refine();
                if (y.is_rational()) return    compare(y.rational());
                if (  is_rational()) return -y.compare(  rational());
                if (  high() <= y.low()) return CGAL::SMALLER;
                if (y.high() <=   low()) return CGAL::LARGER;
            }
        }
#endif
        // filter 2: probabilistically check coprimality
        if (!are_distinct) {
            are_distinct = !(may_have_common_factor(polynomial(),
                                                    y.polynomial()));
        }
        if (!are_distinct) {
            // OK, filters failed. So we have to do the actual work
            // and compute the gcd of the defining polynomials.

            // we have ]low(), high()[ == ]y.low(),y.high()[ == ]L,R[
            // and let both numbers decide for the gcd or its complement
            Poly F1,F2,G;
            G = CGAL::gcd_up_to_constant_factor(polynomial(),y.polynomial()); 
            F1 = CGAL::integral_division_up_to_constant_factor(polynomial(),G);
            CGAL_postcondition(CGAL::degree(F1)==
                               CGAL::degree(polynomial())-CGAL::degree(G));
            F2 = CGAL::integral_division_up_to_constant_factor(y.polynomial(),G);
            CGAL_postcondition(CGAL::degree(F2)==
                               CGAL::degree(y.polynomial())-CGAL::degree(G));
           
            this->learn_from(G,F1);
            y.learn_from(G,F2);

            // this may simplify them due to degree loss
            if (y.is_rational()) return    compare(y.rational());
            if (  is_rational()) return -y.compare(  rational());
        
            // type of x and y is still IS_GENERAL
            // check for equality
            if (G.sign_at(L)!=G.sign_at(R)){
                this->introduce(y);
                return CGAL::EQUAL;
            }
        }
        
        // if we are here, we know the numbers to be distinct
        // refine to disjointness
        for (;;) {
            y.refine();
            refine();
            if (y.is_rational()) return    compare(y.rational());
            if (  is_rational()) return -y.compare(  rational());
            if (  high() <= y.low()) return CGAL::SMALLER;
            if (y.high() <=   low()) return CGAL::LARGER;
        }
    }


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
  
    void refine_poly_approximation() const {
        CGAL_precondition(current_prec > 1);
        current_prec *= 2;
        // std::cout <<"ALGREAL: refine approx: "<<  current_prec<<std::endl;
        set_precision(BFI(),current_prec);
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
    
    void update_poly_approximation() const {
        long old_prec = set_precision( BFI(), current_prec );
        
        polynomial_approx.clear();
        convert_coeffs( this->polynomial(), std::back_inserter( polynomial_approx ) );
                    
        set_precision( BFI(), old_prec );
        
        // TODO: Problems if the next block gets executed
//        Self const *next = static_cast<const Self *>(this->next);
//        while(this != next){
//            // std::cout << this << " " << next << std::endl;
//            next->polynomial_approx = polynomial_approx; 
//            next->current_prec      = current_prec; 
//            next = static_cast<const Self *>(next->next);
//        }
    }
    
    void refine() const{              
        if(this->is_rational()) return;
        
       // std::cout << "refine begin ------- "<< std::endl;

        Field m = (this->low()+this->high())/Field(2);

        // Currently, sign_of_polynomial_at performs exactly the needed
        // refinement. 
        // TODO: But what if this changes? 
        sign_of_polynomial_at( m );
        //std::cout << "refine end ----------------- "<<current_prec<<  std::endl;
    }

protected:
    virtual CGAL::Sign sign_of_polynomial_at( const Field& f ) const {
        //return polynomial().sign_at( f );
        long old_prec = set_precision(BFI(),current_prec);

        Field m = f;
        CGAL::simplify(m);

        // if polynomial has changed, the approximated polynomial gets
        // updated  
        if ( (int) polynomial_approx.size() != 
             CGAL::degree(this->polynomial_)+1) {
            update_poly_approximation();
        }        
        
        CGAL_postcondition(polynomial_approx.size() > 0);
        
        BFI eval = evaluate_polynomial_approx(convert_to_bfi(m));
                
        CGAL::Sign s = CGAL::sign(CGAL::lower(eval));
              
        // correct sign if needed
        if( s*CGAL::sign(CGAL::upper(eval) ) != CGAL::POSITIVE ){
            
            //std::cout << "APPROX FAILED-------------------------------"<<std::endl;
            s = this->polynomial().sign_at(m);
            if ( s != CGAL::ZERO ) {
                refine_poly_approximation(); 
            }
        }
        
        CGAL_postcondition(s == this->polynomial_.sign_at(m));

        if ( s == CGAL::ZERO ) {
            this->learn_from(m);
        }else{
            if ( s == this->sign_at_low() ) this->low_  = m;
            else                            this->high_ = m; 
        }

        set_precision(BFI(),old_prec);
        
        return s;
    }    
    
};

}//namespace internal

} //namespace CGAL

#endif //CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_REP_BFI_H

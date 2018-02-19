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
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_REP_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_REP_H

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Polynomial_traits_d.h>
#include <boost/optional.hpp>

namespace CGAL {

namespace internal {

// definition of the Algebraic_real_rep x:
    
//IS_GENERAL: 
// low_  lower bound of x 
// high_ upper bound of x 
// polynomial_ a square free polynomial 
// sign_at_low_ = polynomial_.evaluate(low_)
// x is the only root of polynomial_ in the open interval ]low_,high_[ 
// low_ != x != high
// ******************* EXEPTION *******************
// x is rational: in this case low=high=x
 
template< class Coefficient_, class Rational_> 
class Algebraic_real_rep{

public:
    typedef Coefficient_                             Coefficient;
    typedef Rational_                                Bound;
    typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type Polynomial_1;

private: 
    typedef Rational_                                Rational;
    typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type Poly;
    typedef CGAL::Sign                              TRI_BOOL;

    typedef Algebraic_real_rep <Coefficient,Rational>     Self;   
public:
    typedef boost::optional< std::pair<double, double> >                  Interval_option; 

    mutable Poly polynomial_;                //!< square free polynomial     
    mutable Rational low_;                      //!< lower endpoint of interval
    mutable Rational high_;                     //!< upper endpoint of interval
    mutable CGAL::Sign sign_at_low_;         //!< sign of polynomial a \a low_
    mutable Interval_option interval_option;               
    mutable bool  is_rational_;

    mutable const Self *next;
    mutable const Self *prev; 
private:
    void copy_all_members(const Algebraic_real_rep& y) const {
        polynomial_     = y.polynomial_;
        low_            = y.low_;
        high_           = y.high_;
        sign_at_low_    = y.sign_at_low_;
        interval_option = y.interval_option;
        is_rational_    = y.is_rational_;
    }
     
    
protected:
    virtual CGAL::Sign sign_of_polynomial_at( const Rational& f ) const {
        return polynomial().sign_at( f );
    }
    
protected:
    // joins the two lists of related algebraic reals  
    void introduce(const Algebraic_real_rep& y) const{
        y.prev->next=next;
        next->prev=y.prev;
        next=&y;
        y.prev=this;
    }
protected:
   
    // x leaves the list of related algebraic reals
    void erase_from_list() const {
        next->prev=prev;
        prev->next=next;
        next=prev=this;
    }
protected:
    virtual void set_explicit_rep(const Rational& m) const {     

        typename Fraction_traits<Rational>::Decompose decomp;
        typename Fraction_traits<Rational>::Numerator_type num;
        typename Fraction_traits<Rational>::Denominator_type den;
        decomp(m, num, den);
        
        polynomial_ = Poly(-Coefficient(num),Coefficient(den)); 
        is_rational_ = true;
        low_  = m;
        high_ = m;
        sign_at_low_= CGAL::ZERO;
        interval_option = Interval_option(CGAL::to_interval(m));
    }
protected:
    virtual void set_implicit_rep(const Poly & P, 
                          const Rational& LOW, 
                          const Rational& HIGH,
                          bool use_expensive_sign = false ) const {
        
        CGAL_precondition(LOW < HIGH);        
        CGAL_precondition(P.sign_at(LOW)  != CGAL::ZERO);       
        CGAL_precondition(P.sign_at(HIGH) != CGAL::ZERO);     
        CGAL_precondition(P.sign_at(HIGH) != P.sign_at(LOW));
        CGAL_precondition_msg(typename CGAL::Polynomial_traits_d< Poly >::Is_square_free()(P), "Polynomial not square-free.");
               
        polynomial_ = P;
        low_=LOW;
        high_=HIGH;
        sign_at_low_=(use_expensive_sign ? P.sign_at(LOW) : sign_of_polynomial_at(LOW) );

        is_rational_=false;

        // interval_option left out

        // trys to set rational if degree is 1 
        typedef typename CGAL::Coercion_traits< Coefficient, Rational >::Type RET;
        set_rational(RET());
    }
private: 
    template < class NTX >
    void set_rational(NTX) const {
    }
    
    void set_rational(Rational) const {
        if (CGAL::degree(polynomial()) == 1) {
            low_ = high_ = -Rational(polynomial()[0])/Rational(polynomial()[1]);
            is_rational_ = true;
            interval_option = Interval_option(CGAL::to_interval(rational()));
        }
    }
    

public:
    //! creates the algebraic real from int \a i.
    explicit Algebraic_real_rep(int i = 0){
        next=prev=this;
        set_explicit_rep(i);
    } 
    //! creates the algebraic real from Rational \a m.
    explicit Algebraic_real_rep(const Rational& m){
        next=prev=this;
        set_explicit_rep(m);        
    } 
    /*! \brief creates the algebraic real as the unique root of \a P
        in the open interval <var>]low,high[</var>.
        \pre the polynomial \a P is square free 
        \pre P(low)!=0
        \pre P(high)<0
        \pre x is the one and only root in the open interval of \a P.
    */
    Algebraic_real_rep(const Poly& P, Rational LOW, Rational HIGH) { 
        if (LOW == HIGH) {
            CGAL_precondition(P.sign_at(LOW) == CGAL::ZERO);
            set_explicit_rep(LOW);        
        } else {
            CGAL_precondition(LOW < HIGH);        
            CGAL_precondition(P.sign_at(LOW)  != CGAL::ZERO);       
            CGAL_precondition(P.sign_at(HIGH) != CGAL::ZERO);     
            CGAL_precondition(P.sign_at(HIGH) != P.sign_at(LOW));
            CGAL_precondition_msg(typename CGAL::Polynomial_traits_d<Poly>::Is_square_free()(P), "Polynomial not square-free.");
            set_implicit_rep(P,LOW,HIGH);
        }
        next=prev=this;
    }

    //! copy constructor
    Algebraic_real_rep(const Self& y){         
        next=prev=this;
        copy_all_members(y);
        introduce(y);  
    }
    //! destructor
    virtual ~Algebraic_real_rep(){
        erase_from_list();
    }

    // assignment 
    Algebraic_real_rep& operator=(const Self& y) {
        if ( this != & y) {
            erase_from_list();
            copy_all_members(y);
            introduce(y);
        }
        CGAL_expensive_postcondition(self_test());
        return *this;
    }

public:
    const Rational&           low()           const {return low_;}
    const Rational&           high()          const {return high_;}
    const Poly&            polynomial()    const {return polynomial_;}
        
    const TRI_BOOL&        sign_at_low()   const {return sign_at_low_;}
    
    /*const std::pair< double, double>& interval() const {
        if( interval_option ) {
            return *interval_option;
        } else {
            typedef typename Get_arithmetic_kernel< Coefficient >::Arithmetic_kernel::Bigfloat BF;
            //typedef typename LEDA_arithmetic_kernel::Bigfloat BF;
            long old_precision = get_precision( BF() );
            set_precision( BF(), 53 );
            interval_option = Interval_option( internal::to_interval( convert_to_bfi( (*this) ) ) );
            set_precision( BF(), old_precision );
            return *interval_option;            
        }
    }*/
    /*const std::pair<double, double>&        interval()      const {
        if(interval_option){
            return *interval_option;
        }else{
            interval_option = Interval_option(
                compute_interval_approximation(polynomial(),low(),high()));
            return *interval_option;
        }
    }*/
    bool                   is_rational()   const {return is_rational_;}
    const Rational&           rational()      const {
        CGAL_precondition(is_rational());
        CGAL_precondition(low()==high());
        CGAL_precondition(polynomial().sign_at(low())==CGAL::ZERO);
        return low_;
    }

    /*! \brief convert algebraic real to type double */
    /*double to_double() const {
        return (to_interval().first +to_interval().second)/2;
    }*/
 
    /*! \brief convert algebraic real to type Interval */
    /*std::pair<double, double> to_interval() const { 
        
        CGAL_precondition(compare(Rational(interval().first)) != CGAL::SMALLER);
        CGAL_precondition(compare(Rational(interval().second)) != CGAL::LARGER);
        return interval();    
    }*/
    
private:
    bool self_test() const{
        if(!is_rational()){
            CGAL_precondition(low()<high());
            
            CGAL_precondition(sign_at_low() != CGAL::ZERO);
            CGAL_precondition(polynomial().sign_at(low() )== sign_at_low());
            CGAL_precondition(polynomial().sign_at(high())==-sign_at_low());
        }else{
            CGAL_precondition(low()==high());             
            CGAL_precondition(polynomial().sign_at(low() )==CGAL::ZERO);         
        }
        return true;
    }    
public:

    virtual void bisect() const {
        
        if(is_rational()) return;
        
        Rational m = (low()+high())/Rational(2);
        CGAL::simplify(m);
        
        CGAL::Sign s = sign_of_polynomial_at(m);
        if ( s == CGAL::ZERO ) learn_from(m);
        else 
            if ( s == sign_at_low() ) low_  = m;
            else                      high_ = m;  
    }

    virtual void refine() const{ 
        bisect();
    }

public:
    template <class NTX>
    void strong_refine(const NTX& m) const{
        
        if(is_rational()) return;

        if( NTX(low()) <= m && m <= NTX(high()) ){
            CGAL_precondition(polynomial().sign_at(m)!=CGAL::ZERO);
            refine();
            while( NTX(low()) <= m && m <= NTX(high())) refine();
        }
    }
public: 
    virtual void refine_at(const Rational& m) const{
        if(is_rational()) return;
        if( m <= low() || high() <= m ) return;
        
        // now: low < m < high 
        CGAL::Sign s = sign_of_polynomial_at(m);
        if(s == CGAL::ZERO)         learn_from(m);
        else (s == sign_at_low())? low_ = m : high_=m;
    }
    
public:     
    CGAL::Comparison_result
    compare(const Rational& y, bool are_distinct = false) const {  
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
                learn_from(y);
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
    compare (const Algebraic_real_rep& y, bool are_distinct = false) const{
        if(  is_rational()) return -y.compare(  rational());
        if(y.is_rational()) return    compare(y.rational());
        
        // type of both x and y IS_GENERAL
        if ( high() <= y.low() ) return CGAL::SMALLER;
        if ( low()  >= y.high()) return CGAL::LARGER;

        // intersection isolating intervals is ]L,R[
        Rational L = (low()  > y.low() ) ? low()  : y.low()  ;
        Rational R = (high() < y.high()) ? high() : y.high() ;

        // refine to smaller intervals at intersection interval boundaries
        // this can change type() only to IS_RATIONAL
        refine_at(L);
        refine_at(R);
        y.refine_at(L);
        y.refine_at(R);

        if (  is_rational()) return -y.compare(  rational());
        if (y.is_rational()) return    compare(y.rational());

        // type of both x and y still IS_GENERAL
        if ( high()   <= y.low() ) return CGAL::SMALLER;
        if ( y.high() <=   low() ) return CGAL::LARGER;
       
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
            F1 = integral_division_up_to_constant_factor(polynomial(),G);
            CGAL_postcondition(CGAL::degree(F1)==
                               CGAL::degree(polynomial())-CGAL::degree(G));
            F2 = integral_division_up_to_constant_factor(y.polynomial(),G);
            CGAL_postcondition(CGAL::degree(F2)==
                               CGAL::degree(y.polynomial())-CGAL::degree(G));
           
            learn_from(G,F1);
            y.learn_from(G,F2);

            // this may simplify them due to degree loss
            if (y.is_rational()) return    compare(y.rational());
            if (  is_rational()) return -y.compare(  rational());
        
            // type of x and y is still IS_GENERAL
            // check for equality
            if (G.sign_at(L)!=G.sign_at(R)){
                introduce(y);
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

protected:
    void learn_from(const Rational& m) const {
        if(is_rational()){
            CGAL_precondition(rational()==m);
            return;
        }
        CGAL_precondition(polynomial().sign_at(m)==CGAL::ZERO);
        typename Fraction_traits<Rational>::Decompose decomp;
        typename Fraction_traits<Rational>::Numerator_type num;
        typename Fraction_traits<Rational>::Denominator_type den;
        decomp(m, num, den);
        Poly G = Poly(-Coefficient(num),Coefficient(den)); 
        Poly F1= integral_division_up_to_constant_factor(polynomial(),G);
        CGAL_postcondition(CGAL::degree(F1)==
                           CGAL::degree(polynomial())-CGAL::degree(G));
   
        learn_from(G,F1);
     
        if(!is_rational()) set_explicit_rep(m);
    }   
  
    // splits the list of related algebraic reals according to the two
    // factors of the polynomial
    void learn_from(const Poly& pfactor1,
                    const Poly& pfactor2) const {
        if(CGAL::degree(pfactor1)==0) return;
        if(CGAL::degree(pfactor2)==0) return;
         
        CGAL_precondition(may_have_common_factor(polynomial(),pfactor1));
        CGAL_precondition(may_have_common_factor(polynomial(),pfactor2));
        CGAL_precondition(CGAL::degree(pfactor1)+CGAL::degree(pfactor2)
                          == CGAL::degree(polynomial()));


      
        Self dummy1, dummy2;
        Self const *last=prev;
        Self const *current=this;
        Self const *next;
        while(true){
            next=current->next;
            
            if(pfactor1.sign_at(current->low())!=
               pfactor1.sign_at(current->high())){
                current->next=&dummy1;
                current->prev=dummy1.prev;
                dummy1.prev->next=current;
                dummy1.prev=current;
                current->set_implicit_rep(pfactor1,
                                          current->low(),
                                          current->high(), true );
            }else{
                current->next=&dummy2;
                current->prev=dummy2.prev;
                dummy2.prev->next=current;
                dummy2.prev=current;                    
                current->set_implicit_rep(pfactor2,
                                          current->low(),
                                          current->high(), true );
            }
            
            if(current==last) break;
            current=next;
        }
    }
    
public:
    bool is_root_of(const Poly& Q) const {
        if (CGAL::degree(Q) == 0) return Q.is_zero();
        if (is_rational()  ) return CGAL::ZERO == Q.sign_at(rational());
        
        if ( may_have_common_factor(polynomial(), Q) ) {
          Poly G = CGAL::gcd_up_to_constant_factor(polynomial(),Q);
            if(CGAL::degree(G)!=0){
                Poly F1 = integral_division_up_to_constant_factor(
                        polynomial(),G
                );
                CGAL_postcondition(CGAL::degree(F1)==
                                   CGAL::degree(polynomial())-CGAL::degree(G));
                learn_from(G,F1);
                return (G.sign_at(low()) != G.sign_at(high()));
            }
        } 
        return false;      
    }
public: 
    virtual void simplify() const{
        if(is_rational()){
          
        }else{
            CGAL::simplify(low_);
            CGAL::simplify(high_);
            polynomial_ = typename Polynomial_traits_d<Poly>::Canonicalize()(polynomial_);
        }
    }
};
} // namespace internal

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_REP_H

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
//                  Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_ALGEBRAIC_REAL_QUADRATIC_REFINEMENT_REP_BFI_H
#define CGAL_ALGEBRAIC_REAL_QUADRATIC_REFINEMENT_REP_BFI_H

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep.h>

#include <CGAL/Algebraic_kernel_d/Float_traits.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/convert_to_bfi.h>

#include <CGAL/Arithmetic_kernel.h>

#include <boost/none.hpp>

#include <CGAL/Polynomial_type_generator.h>

namespace CGAL {

namespace internal {

// definition of the Algebraic_real_rep x:
    
//For details about the method, see
/*
 * @Unpublished{abbott-quadratic,
 *    author = 	 {John Abbott},
 *    title = 	 {Quadratic Interval Refinement for Real Roots},
 *    url =      {http://www.dima.unige.it/~abbott/},
 *  note =       {Poster presented at the 2006 Internat. Sympos. on Symbolic
 and Algebraic Computation (ISSAC 2006)}
 * }
 */

template< class Coefficient_, class Field_> 
class Algebraic_real_quadratic_refinement_rep_bfi
    : public Algebraic_real_rep<Coefficient_, Field_> {

    typedef Coefficient_                            Coefficient;
    typedef Field_                                  Field;

    typedef typename 
    CGAL::Get_arithmetic_kernel<Field>::Arithmetic_kernel:: 
        Bigfloat_interval BFI;

    typedef typename CGAL::Bigfloat_interval_traits<BFI>::Bound BF;

    // This is a implicit restriction - Field must be some type
    // modelling rational numbers to get an integer type
    typedef typename CGAL::Fraction_traits<Field>::Numerator_type Integer;
    
    typedef typename
        CGAL::Polynomial_type_generator<Coefficient,1>::Type Poly;

    typedef Algebraic_real_rep <Coefficient,Field>     Base;
    typedef Algebraic_real_quadratic_refinement_rep_bfi<Coefficient,Field> 
    Self;

    typedef typename CGAL::Coercion_traits<Coefficient,Field>::Type 
    Eval_result_type;

private:

    mutable long prec_;

    typedef typename 
        CGAL::Polynomial_traits_d<Poly>::template Rebind<BFI,1>
    ::Other::Type BFI_polynomial;

    mutable boost::optional
        < BFI_polynomial > f_bfi_;
    
    mutable boost::optional<BFI> low_bfi_, f_low_bfi_, 
        high_bfi_, f_high_bfi_; 

    mutable long N;

  // TODO: replace by call of Coercion_traits<Poly,BFI>::Cast() 
    BFI_polynomial _convert_polynomial_to_bfi(const Poly& f) const {
        std::vector<BFI> coeffs;
        for(int i = 0; i <= CGAL::degree(f); i++) {
            coeffs.push_back(CGAL::convert_to_bfi(f[i]));
        }
        return BFI_polynomial(coeffs.begin(), coeffs.end());   
    }

    void _set_prec(long new_prec) const {
        
        prec_ = new_prec;
        CGAL::set_precision(BFI(), prec_);

        f_bfi_ = _convert_polynomial_to_bfi(this->polynomial());

        low_bfi_ = CGAL::convert_to_bfi(this->low());
        
        high_bfi_ = CGAL::convert_to_bfi(this->high());
        f_low_bfi_ = f_bfi_.get().evaluate(low_bfi_.get());
        f_high_bfi_ = f_bfi_.get().evaluate(high_bfi_.get());

    }

    CGAL::Sign _sign_at(Field m, BFI& m_bfi, BFI& f_m_bfi) const {

        if(! f_bfi_) {
            f_bfi_ = _convert_polynomial_to_bfi(this->polynomial());
        }

        m_bfi = CGAL::convert_to_bfi(m);
        f_m_bfi = f_bfi_.get().evaluate(m_bfi);
        
        if(CGAL::zero_in(f_m_bfi)) {
            
            // Okay, compute exactly
            return CGAL::sign(this->polynomial().evaluate(m));
        }
        
        // If we are here, then the interval is away from zero
        
        return CGAL::sign(CGAL::upper(f_m_bfi));

    }


private:
    // Stores whether the last bisection has taken the upper or lower part
    mutable bool last_bisect_lower;

public:
    //! creates the algebraic real from int \a i.
    explicit Algebraic_real_quadratic_refinement_rep_bfi(int i = 0)
        : Base(i),N(2){
    } 
    //! creates the algebraic real from Field \a m.
    explicit Algebraic_real_quadratic_refinement_rep_bfi(const Field& m)
        : Base(m),N(2) {
    } 
    /*! \brief creates the algebraic real as the unique root of \a P
      in the open interval <var>]low,high[</var>.
      \pre the polynomial \a P is square free 
      \pre P(low)!=0
      \pre P(high)<0
      \pre x is the one and only root in the open interval of \a P.
    */
    Algebraic_real_quadratic_refinement_rep_bfi(const Poly& P, 
                                                Field LOW, 
                                                Field HIGH) 
        : Base(P,LOW,HIGH), 
          N(2)
    { 
        _set_prec(16);
        
    }

    //! copy constructor
    Algebraic_real_quadratic_refinement_rep_bfi(const Self& y)
        : Base(y), prec_(y.prec_), f_bfi_(y.f_bfi_), low_bfi_(y.low_bfi_),
          f_low_bfi_(y.f_low_bfi_), high_bfi_(y.high_bfi_),
          f_high_bfi_(y.f_high_bfi_), N(y.N)
    {
    }

    // assignment 
    Algebraic_real_quadratic_refinement_rep_bfi& operator=(const Self& y) {
        Base::operator=(y);
	f_low_bfi_=y.f_low_bfi_;
	f_high_bfi_=y.f_high_bfi_;
        low_bfi_=y.low_bfi_;
	high_bfi_=y.high_bfi_;
        f_bfi_=y.f_bfi_;
	N=y.N;
        return *this;
    }

public:

    virtual void bisect() const{

        if(this->is_rational()) return;
        
        Field m = (this->low_+this->high_)/Field(2);

        CGAL::simplify(m);
        BFI m_bfi, f_m_bfi;
        CGAL::Sign s = _sign_at(m,m_bfi,f_m_bfi);

        if (s  == ::CGAL::ZERO ) {
            this->learn_from(m);
	}
        else {
            if ( s == this->sign_at_low() ) {
                this->low_  = m;
                low_bfi_ = m_bfi;
                f_low_bfi_ = f_m_bfi;
                last_bisect_lower=false;
            } else {
                this->high_ = m; 
                high_bfi_ = m_bfi;
                f_high_bfi_ = f_m_bfi;
                last_bisect_lower=true;
            }
        }

    }

protected:
    virtual void set_implicit_rep(const Poly & P, 
				  const Field& LOW, 
				  const Field& HIGH,
                                  bool dummy_bool=false) const {

        bool poly_changed = (P!=this->polynomial());
        if(poly_changed) {
            f_bfi_ = boost::none;
        }
        if(poly_changed || LOW != this->low()) {
            f_low_bfi_ = low_bfi_ = boost::none;
        }
        if(poly_changed || HIGH != this->high()) {
            f_high_bfi_ = high_bfi_ = boost::none;
        }
        Base::set_implicit_rep(P,LOW,HIGH,dummy_bool);
    }

    virtual void set_explicit_rep(const Field& m) const {
        f_bfi_ = boost::none;
        f_low_bfi_ = low_bfi_ = boost::none;
        f_high_bfi_ = high_bfi_ = boost::none;
        Base::set_explicit_rep(m);
    }

     
public: 
    virtual void refine_at(const Field& m) const{
        Field old_low_=this->low_, old_high_=this->high_;
        Poly old_pol = this->polynomial();
        Base::refine_at(m);
        if(this->is_rational()) return;

        if(old_low_!=this->low_) {
            f_low_bfi_ = low_bfi_ = boost::none;
        }
        if(old_high_!=this->high_) {
            f_high_bfi_ = high_bfi_ = boost::none;
        }
        if(old_pol != this->polynomial()) {
            f_bfi_ = boost::none;
        }
    }
    
    // Abbott's refinement method
    virtual void refine() const {

        if(this->is_rational()) {
            return;
        }

        CGAL_assertion(this->low() != this->high());

        long old_prec = CGAL::get_precision(BFI());

        CGAL::set_precision(BFI(),prec_);

        CGAL_assertion(CGAL::sign(this->polynomial().evaluate(this->low()))
                       ==this->sign_at_low_);

        CGAL_assertion( this->sign_at_low_ != 
                        CGAL::sign(this->polynomial().evaluate(this->high())) );
      
        CGAL_assertion(this->low() != this->high());
      
        Integer i = find_interval();
      
        while(N!=1 && !refine_by_factor(i)) {
            N/=2;
            i = find_interval();
        }
        N*=2;

        CGAL::set_precision(BFI(),old_prec);

    }

private:


    std::pair<Integer,Integer> _to_integer_interval(BFI z, long N) const {

        Integer i_low, i_high;

        //typename CGAL::internal::Real_embeddable_extension<BF>::Floor floor;
        //typename CGAL::internal::Real_embeddable_extension<BF>::Ceil ceil;
        typename CGAL::internal::Float_traits<BF>::Mul_by_pow_of_2 mul_2;

        BF z_low=CGAL::lower(z), z_high = CGAL::upper(z);

        i_low = CGAL::internal::floor(mul_2(CGAL::lower(z),N));
        i_high = CGAL::internal::ceil(mul_2(CGAL::upper(z),N));

        return std::make_pair(i_low,i_high);

    }

    Integer find_interval() const {

        if(! f_bfi_) {
            f_bfi_ = _convert_polynomial_to_bfi(this->polynomial());
        }

        if(! low_bfi_) {
            low_bfi_ = CGAL::convert_to_bfi(this->low());
        }
        if(! f_low_bfi_) {
            f_low_bfi_ = f_bfi_.get().evaluate(low_bfi_.get());
        }
        if(! high_bfi_) {
            high_bfi_ = CGAL::convert_to_bfi(this->high());
        }
        if(! f_high_bfi_) {
            f_high_bfi_ = f_bfi_.get().evaluate(high_bfi_.get());
        }
        Integer i;
        while(true) {

            if(CGAL::zero_in(f_low_bfi_.get() - f_high_bfi_.get())) {
                _set_prec(2*prec_);
                continue;
            }

            BFI denom = f_low_bfi_.get()-f_high_bfi_.get();

            BFI z = f_low_bfi_.get() / denom;

            std::pair<Integer, Integer> int_pair = _to_integer_interval(z,N);
            Integer i_low = int_pair.first;
            Integer i_high = int_pair.second;

            if(CGAL::abs(i_high-i_low) <= 2) {
                i = CGAL::abs((i_high+i_low))/2;
                break;
            }
            _set_prec(2*prec_);
                        
        }

        CGAL_postcondition(i>=0 && 
                           i <= CGAL::ipower(Integer(2),N));
        
        return i;
    }

    bool refine_by_factor(Integer i) const {

        if(N==2) {
            bool refined = refine_by_factor_4(i);
            return refined;
        }
        bool refined = refine_by_factor_greater_4(i);
        return refined;
    }

    bool refine_by_factor_4(Integer i) const {
        Integer actual_i;
        bisect();
        actual_i=last_bisect_lower ? 0 : 2;
        bisect();
        if(! last_bisect_lower) {
            actual_i = actual_i+1;
        }
        return this->is_rational() || actual_i==i;
    }

    bool refine_by_factor_greater_4(Integer i) const {
        Integer intervals = CGAL::ipower(Integer(2),N);
        Field step = (this->high_-this->low_)/Field(intervals);
        CGAL::simplify(step);
        Field m
            = this->low_ + step*Field(i);
        CGAL::simplify(m);
        BFI m_bfi, f_m_bfi;
        CGAL::Sign s_m = _sign_at(m,m_bfi,f_m_bfi);
        if(s_m==CGAL::ZERO) {
            this->learn_from(m);
            return true;
        }
        Field new_left, new_right;
        BFI new_left_bfi, new_right_bfi, f_new_left_bfi, f_new_right_bfi;
        CGAL::Sign s_new_left, s_new_right;
        if(s_m == this->sign_at_low_) {
            // Go to the right
            new_left=m;
            new_left_bfi = m_bfi;
            f_new_left_bfi = f_m_bfi;
            s_new_left = s_m;
        
            new_right = m+step;
            CGAL::simplify(new_right);
            s_new_right = _sign_at(new_right,new_right_bfi,f_new_right_bfi);
            if(s_new_right==CGAL::ZERO) {
                this->learn_from(new_right);
                return true;
            }
        } else {
            // Go to the left
            new_right=m;
            new_right_bfi = m_bfi;
            f_new_right_bfi = f_m_bfi;
            s_new_right=s_m;
            new_left = m-step;
            CGAL::simplify(new_left);
            s_new_left = _sign_at(new_left,new_left_bfi,f_new_left_bfi);
            if(s_new_left==CGAL::ZERO) {
                this->learn_from(new_left);
                return true;
            }
        }
        if(s_new_left != s_new_right) {
            this->low_=new_left;
            this->high_=new_right;
            low_bfi_ = new_left_bfi;
            high_bfi_ = new_right_bfi;
            f_low_bfi_ = f_new_left_bfi;
            f_high_bfi_ = f_new_right_bfi;
            this->sign_at_low_=s_new_left;
            return true;
        }
        else {
            return false;
        }
    }

protected:
    virtual CGAL::Sign sign_of_polynomial_at( const Field& f ) const {
        //return polynomial().sign_at( f );

        Field m = f;
        CGAL::simplify(m);

        if(! f_bfi_) {
            f_bfi_ = _convert_polynomial_to_bfi(this->polynomial());
        }

        BFI eval = f_bfi_.get().evaluate(convert_to_bfi(m));
        
        CGAL::Sign s = CGAL::sign(CGAL::lower(eval));
              
        // correct sign if needed
        if( s*CGAL::sign(CGAL::upper(eval) ) != CGAL::POSITIVE ){
            
            //std::cout << "APPROX FAILED-------------------------------"<<std::endl;
            s = this->polynomial().sign_at(m);
            if ( s != CGAL::ZERO ) {
                _set_prec(2*prec_);
            }
        }
        

        CGAL_postcondition(s == this->polynomial_.sign_at(m));

        if ( s == CGAL::ZERO ) {
            this->learn_from(m);
        }else{
            if ( s == this->sign_at_low() ) this->low_  = m;
            else                            this->high_ = m; 
        }

        return s;
    }        

public: 
    virtual void simplify() const {
        Poly f_old = this->polynomial();
        Base::simplify();
        if(f_old != this->polynomial()) {
            f_bfi_ = boost::none;
        }      
    }
};
} // namepace internal

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_ALGEBRAIC_REAL_QUADRATIC_REFINEMENT_REP_H

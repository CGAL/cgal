// Copyright (c) 2006-2008 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_ALGEBRAIC_1_MEMBER_H
#define CGAL_RS_ALGEBRAIC_1_MEMBER_H

#include <CGAL/RS/refine_1_rs.h>
#include <CGAL/assertions.h>

namespace CGAL{

inline
mpfi_srcptr Algebraic_1::mpfi()const{
        return Ptr()->_mpfi;
}

inline
mpfi_ptr Algebraic_1::mpfi(){
        return ptr()->_mpfi;
}

inline
Gmpfi Algebraic_1::interval()const{
        return Gmpfi(Ptr()->_mpfi);
}

inline
Gmpfr Algebraic_1::inf()const{
        return Gmpfr(&(mpfi()->left));
}

inline
Gmpfr Algebraic_1::sup()const{
        return Gmpfr(&(mpfi()->right));
}

inline
const RS_polynomial_1& Algebraic_1::pol()const{
        return *(Ptr()->_poly);
}

inline
int Algebraic_1::nr()const{
        return ptr()->_nr;
}

inline
int Algebraic_1::mult()const{
        return ptr()->_mult;
}

inline
void Algebraic_1::set_mpfi_ptr(mpfi_srcptr x){
        // *mpfi()=*x;
        // mpfi_set(mpfi(),x);
        set_mpfi(x);
}

inline
void Algebraic_1::clear_pol(){
        ptr()->_poly=NULL;
}

inline
void Algebraic_1::set_pol(const RS_polynomial_1 &p){
        ptr()->_poly=const_cast<RS_polynomial_1*>(&p);
}

inline
void Algebraic_1::set_nr(int n){
        ptr()->_nr=n;
}

inline
void Algebraic_1::set_mult(int m){
        ptr()->_mult=m;
}

inline
void Algebraic_1::set_prec(mp_prec_t p){
        mpfi_round_prec(mpfi(),p);
}

inline
void Algebraic_1::set_lefteval(Sign s)const{
        Ptr()->_lefteval=s;
}

inline
mp_prec_t Algebraic_1::get_prec()const{
        return mpfi_get_prec(mpfi());
}

inline
mpfr_srcptr Algebraic_1::left()const{
        return &(mpfi()->left);
}

inline
mpfr_srcptr Algebraic_1::right()const{
        return &(mpfi()->right);
}

inline
Sign Algebraic_1::lefteval()const{
        return ptr()->_lefteval;
}

inline
bool Algebraic_1::is_consistent()const{
        return(&pol()==NULL?false:true);
}

inline
bool Algebraic_1::is_point()const{
        return(mpfr_equal_p(&(mpfi()->left),&(mpfi()->right))!=0);
}

inline
bool Algebraic_1::contains(int n)const{
        return(mpfr_cmp_si(&(mpfi()->left),n)<=0 &&
                        mpfr_cmp_si(&(mpfi()->right),n)>=0);
}

inline
bool Algebraic_1::contains(mpfr_srcptr n)const{
        return(mpfr_lessequal_p(&(mpfi()->left),n) &&
                        mpfr_greaterequal_p(&(mpfi()->right),n));
}

inline
bool Algebraic_1::contains(const Gmpz &n)const{
        return(mpfr_cmp_z(&(mpfi()->left),n.mpz())<=0 &&
                        mpfr_cmp_z(&(mpfi()->right),n.mpz())>=0);
}

inline
std::pair<double,double> Algebraic_1::to_interval()const{
        return std::make_pair(
                        mpfr_get_d(left(),GMP_RNDD),
                        mpfr_get_d(right(),GMP_RNDU));
}

inline
void Algebraic_1::set_mpfi(mpfi_srcptr x){
        mp_prec_t xp;
        xp=mpfi_get_prec(x);
        if(xp>mpfr_get_prec(left()) || xp>mpfr_get_prec(right()))
                mpfi_set_prec(mpfi(),xp);
        mpfi_set(mpfi(),x);
}

inline
bool Algebraic_1::overlaps(const Algebraic_1&a)const{
        if(mpfr_lessequal_p(left(),a.left()))
                return (mpfr_lessequal_p(a.left(),right())!=0);
        else
                return (mpfr_lessequal_p(left(),a.right())!=0);
}

inline
bool Algebraic_1::is_valid()const{
        return (mpfi_nan_p(mpfi())==0);
}

inline
bool Algebraic_1::is_finite()const{
        return (mpfi_bounded_p(mpfi())!=0);
}

//template <class _Gcd_policy>
inline
double Algebraic_1::to_double()const{
        /*typedef _Gcd_policy     Gcd;
        while(mpfr_get_d(left(),GMP_RNDU)!=mpfr_get_d(right(),GMP_RNDU))
                bisect_n<Gcd>(*this,33);*/
        RS3::refine_1(*this,100);
        CGAL_assertion(mpfr_get_d(left(),GMP_RNDD)==
                       mpfr_get_d(right(),GMP_RNDD));
        CGAL_assertion(mpfr_get_d(left(),GMP_RNDU)==
                       mpfr_get_d(right(),GMP_RNDU));
        CGAL_assertion(mpfr_get_d(left(),GMP_RNDN)==
                       mpfr_get_d(right(),GMP_RNDN));
        return mpfr_get_d(right(),GMP_RNDU);
}

inline
Algebraic_1 Algebraic_1::sqrt()const{
        mpfi_t s;
        mpfi_init(s);
        mpfi_sqrt(s,mpfi());
        Algebraic_1 ret(s);
        return ret;
}

} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_1_MEMBER_H

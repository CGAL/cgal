// Copyright (c) 2006-2009 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_POLYNOMIAL_1_MEMBER_H
#define CGAL_RS_POLYNOMIAL_1_MEMBER_H

namespace CGAL{

//////////////////////
// private functions

inline
void RS_polynomial_1::create_storage(int size){
        _coef=(mpz_t*)(*_allocf)(sizeof(mpz_t)*size);
        _capacity=size;
}

inline
void RS_polynomial_1::free_storage(){
        void *(*af)(size_t);
        void *(*rf)(void*,size_t,size_t);
        void (*ff)(void*,size_t);
        mp_get_memory_functions(&af,&rf,&ff);
        mp_set_memory_functions(_allocf,_reallocf,_freef);
        for(int i=0;i<_degree+1;++i)
                mpz_clear(coef(i));
        mp_set_memory_functions(af,rf,ff);
        (*_freef)(_coef,sizeof(mpz_t)*_capacity);
        _capacity=0;
}

inline
void RS_polynomial_1::fetch_gmp_functions(){
        mp_get_memory_functions(&_allocf,&_reallocf,&_freef);
}

//////////////////////
// private functions

// dangerous function! this will erase the coefficients if the new coefficients
// don't fit in the previously allocated space; otherwise, it will set to zero
// the new coefficients
inline
void RS_polynomial_1::set_degree(int d){
        if(d+1>_capacity){
                free_storage();
                _degree=d;
                ++d;
                create_storage(d);
                for(int i=0;i<=_degree;++i)
                        mpz_init(coef(i));
        }else{
                _degree=d;
                for(int i=_degree+1;i<=d;++i)
                        mpz_init(coef(i));
        }
        return;
}

inline
void RS_polynomial_1::force_degree(int n){
        _degree=n;
}

// to change the storage capacity of a polynomial object
inline
int RS_polynomial_1::resize(int newcap){
        if(newcap<=_capacity)
                return -1;
        int i;
        mpz_t *newcoef=(mpz_t*)(*_allocf)(newcap*sizeof(mpz_t));
        for(i=0;i<=_degree;++i){
                mpz_init_set(newcoef[i],_coef[i]);
                mpz_clear(_coef[i]);
        }
        for(i=_degree+1;i<newcap;++i)
                mpz_init(newcoef[i]);
        (*_freef)(_coef,sizeof(mpz_t)*_capacity);
        _coef=newcoef;
        _capacity=newcap;
        return newcap;
}

inline
void RS_polynomial_1::set_coef(int pow_x,mpz_srcptr z){
        mpz_set(_coef[pow_x],z);
}

inline
void RS_polynomial_1::set_coef(int pow_x,const CGAL::Gmpz &z){
        mpz_set(_coef[pow_x],z.mpz());
}

inline
void RS_polynomial_1::set_coef_ui(int pow_x,unsigned long z){
        mpz_set_ui(_coef[pow_x],z);
}

inline
void RS_polynomial_1::set_coef_si(int pow_x,long z){
        mpz_set_si(_coef[pow_x],z);
}

inline
int RS_polynomial_1::get_degree()const{
        while(!mpz_sgn(coef(_degree))&&_degree)
                --_degree;
        return _degree;
}

inline
int RS_polynomial_1::get_degree_static()const{
        return _degree;
}

inline
bool RS_polynomial_1::has_sfpart()const{
        return (_is_sf?true:(_sfpart.get()!=NULL?true:false));
}

inline
const RS_polynomial_1& RS_polynomial_1::sfpart()const{
        if(_is_sf)
                return *this;
        else
                return *_sfpart;
}

inline
void RS_polynomial_1::set_sfpart(RS_polynomial_1 *s)const{
        _is_sf=false;
        _sfpart=polyptr(s);
}

inline
void RS_polynomial_1::set_sfpart(const polyptr &s)const{
        _is_sf=false;
        _sfpart=s;
}

inline
void RS_polynomial_1::set_sf()const{
        _is_sf=true;
}

inline
bool RS_polynomial_1::has_sqfr()const{
        return (_sqfr.get()!=NULL?true:false);
}

inline
sqfrvec& RS_polynomial_1::sqfr()const{
        return *_sqfr;
}

inline
void RS_polynomial_1::set_sqfr(sqfrvec *s)const{
        _sqfr=sqfrptr(s);
}

inline
void RS_polynomial_1::set_sqfr(const sqfrptr &s)const{
        _sqfr=s;
}

inline
mpz_ptr RS_polynomial_1::leading_coefficient()const{
        return(_coef[get_degree()]);
}

// gets the power of the lowest coefficient non-zero monomial
inline
int RS_polynomial_1::first_non_zero()const{
        int i=0;
        while(i<=_degree)
                if(mpz_sgn(coef(i)))
                        return i;
                else
                        ++i;
        return -1;      // if all the coefficients are zero
}

inline
mpz_t* RS_polynomial_1::get_coefs()const{
        return _coef;
}

inline
mpz_ptr RS_polynomial_1::coef(int i)const{
        return _coef[i];
}

inline
RS_polynomial_1& RS_polynomial_1::derive()const{
        int d=get_degree();
        RS_polynomial_1 *derivative=new RS_polynomial_1(d-1);
        for(int x=0;x<d;++x)
                mpz_mul_si(derivative->coef(x),coef(x+1),x+1);
        return *derivative;
}

inline
RS_polynomial_1& RS_polynomial_1::minusx()const{
        int d=get_degree();
        RS_polynomial_1 *mx=new RS_polynomial_1(d);
        for(int x=0;x<=d;++x){
                if(x%2==1)
                        mpz_set(mx->coef(x),coef(x));
                else
                        mpz_neg(mx->coef(x),coef(x));
        }
        return *mx;
}

// This function copies into lb a number which is less than all roots. It just
// looks for the coefficient with the greatest absolute value. This is useful
// when evaluating curves in the infinite. We can refine this, but it makes no
// sense since it will be slower.
inline
void RS_polynomial_1::get_lower_bound(mpfr_ptr lb)const{
        int d=get_degree();
        mpz_t temp;
        mpz_init(temp);
        mpz_abs(temp,leading_coefficient());
        for(int i=0;i<d;++i)
                if(mpz_cmpabs(coef(i),temp)>0)
                        mpz_abs(temp,coef(i));
        mpz_add_ui(temp,temp,1);
        mpz_neg(temp,temp);
        mpfr_set_z(lb,temp,GMP_RNDD);
        mpz_clear(temp);
}

// the same that above, but the copied value is the upper bound of the roots
inline
void RS_polynomial_1::get_upper_bound(mpfr_ptr ub)const{
        int d=get_degree();
        mpz_t temp;
        mpz_init(temp);
        mpz_abs(temp,leading_coefficient());
        for(int i=0;i<d;++i)
                if(mpz_cmpabs(coef(i),temp)>0)
                        mpz_abs(temp,coef(i));
        mpz_add_ui(temp,temp,1);
        mpfr_set_z(ub,temp,GMP_RNDU);
        mpz_clear(temp);
}

// returns this polynomial times c*x^p
inline
RS_polynomial_1
RS_polynomial_1::times_monomial(mpz_srcptr c,int p)const{
        int dr=get_degree()+p;
        RS_polynomial_1 r(dr);
        for(int i=p;i<=dr;++i)
                mpz_mul(r.coef(i),coef(i-p),c);
        return r;
}

} // namespace CGAL

#endif  // CGAL_RS_POLYNOMIAL_1_MEMBER_H

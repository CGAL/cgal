// Copyright (c) 2007-2009 Inria Lorraine (France). All rights reserved.
// 
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author: Luis Peñaranda <luis.penaranda@loria.fr>

// TODO:
//      - add constructor from (integer,exponent)

#ifndef CGAL_GMPFR_TYPE_H
#define CGAL_GMPFR_TYPE_H

#include <gmp.h>
#include <mpfr.h>
#include <boost/operators.hpp>
#include <CGAL/Handle_for.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <string>

namespace CGAL{

class Gmpfr;

bool operator<(const Gmpfr&,const Gmpfr&);
bool operator==(const Gmpfr&,const Gmpfr&);

bool operator<(const Gmpfr&,long);
bool operator>(const Gmpfr&,long);
bool operator==(const Gmpfr&,long);

bool operator<(const Gmpfr&,unsigned long);
bool operator>(const Gmpfr&,unsigned long);
bool operator==(const Gmpfr&,unsigned long);

bool operator<(const Gmpfr&,int);
bool operator>(const Gmpfr&,int);
bool operator==(const Gmpfr&,int);

bool operator<(const Gmpfr&,double);
bool operator>(const Gmpfr&,double);
bool operator==(const Gmpfr&,double);

bool operator<(const Gmpfr&,long double);
bool operator>(const Gmpfr&,long double);
bool operator==(const Gmpfr&,long double);

bool operator<(const Gmpfr&,const Gmpz&);
bool operator>(const Gmpfr&,const Gmpz&);
bool operator==(const Gmpfr&,const Gmpz&);

bool operator<(const Gmpfr&,const Gmpq&);
bool operator>(const Gmpfr&,const Gmpq&);
bool operator==(const Gmpfr&,const Gmpq&);

struct Gmpfr_rep{
        mpfr_t floating_point_number;
        bool clear_on_destruction;
        Gmpfr_rep():clear_on_destruction(true){}
        ~Gmpfr_rep(){
                if(clear_on_destruction)
                        mpfr_clear(floating_point_number);
        }
};

class Gmpfr:
#ifdef CGAL_GMPFR_NO_REFCOUNT
        Gmpfr_rep,
#else
        Handle_for<Gmpfr_rep>,
#endif
        boost::ordered_euclidian_ring_operators1<Gmpfr,
        boost::ordered_euclidian_ring_operators2<Gmpfr,long,
        boost::ordered_euclidian_ring_operators2<Gmpfr,unsigned long,
        boost::ordered_euclidian_ring_operators2<Gmpfr,int,
        boost::totally_ordered2<Gmpfr,double,
        boost::totally_ordered2<Gmpfr,long double,
        boost::ordered_euclidian_ring_operators2<Gmpfr,Gmpz,
        boost::ordered_euclidian_ring_operators2<Gmpfr,Gmpq
        > > > > > > > >
{
        typedef Handle_for<Gmpfr_rep>   Base;

        public:

        // access

        inline
        mpfr_srcptr fr()const{
#ifdef CGAL_GMPFR_NO_REFCOUNT
                return floating_point_number;
#else
                return Ptr()->floating_point_number;
#endif
        }

        inline
        mpfr_ptr fr(){
#ifdef CGAL_GMPFR_NO_REFCOUNT
                return floating_point_number;
#else
                return ptr()->floating_point_number;
#endif
        }

        inline
        void dont_clear_on_destruction(){
#ifdef CGAL_GMPFR_NO_REFCOUNT
                clear_on_destruction=false;
#else
                ptr()->clear_on_destruction=false;
#endif
        }

        // construction

        Gmpfr(){
                mpfr_init(fr());
        }

        Gmpfr(mpfr_srcptr f){
                mpfr_custom_init_set(
                        fr(),
                        mpfr_custom_get_kind(f),
                        mpfr_custom_get_exp(f),
                        mpfr_get_prec(f),
                        mpfr_custom_get_mantissa(f));
                dont_clear_on_destruction();
                CGAL_assertion(mpfr_equal_p(f,fr()));
        }

        Gmpfr(mpfr_srcptr f,
                        mp_rnd_t r,
                        mp_prec_t p=0){
                mpfr_init2(fr(),p?p:mpfr_get_prec(f));
                mpfr_set(fr(),f,r);
        }

        Gmpfr(mpfr_srcptr f,mp_prec_t p){
                mpfr_init2(fr(),p);
                mpfr_set(fr(),f,get_defrnd());
        }

#define _GMPFR_CONSTRUCTOR_FROM_TYPE(_type,_fun) \
        Gmpfr(_type x,mp_rnd_t r=get_defrnd(),mp_prec_t p=get_defprec()){ \
                mpfr_init2(fr(),p); \
                _fun(fr(),x,r); \
        } \
        Gmpfr(_type x,mp_prec_t p){ \
                mpfr_init2(fr(),p); \
                _fun(fr(),x,get_defrnd()); \
        }

#define _GMPFR_CONSTRUCTOR_FROM_OBJECT(_class,_member,_fun) \
        Gmpfr(          const _class &x, \
                        mp_rnd_t r=get_defrnd(), \
                        mp_prec_t p=get_defprec()){ \
                mpfr_init2(fr(),p); \
                _fun(fr(),x._member,r); \
        } \
        Gmpfr(const _class &z,mp_prec_t p){ \
                mpfr_init2(fr(),p); \
                _fun(fr(),z._member,get_defrnd()); \
        }

        _GMPFR_CONSTRUCTOR_FROM_TYPE(int,mpfr_set_si);
        _GMPFR_CONSTRUCTOR_FROM_TYPE(long,mpfr_set_si);
        _GMPFR_CONSTRUCTOR_FROM_TYPE(unsigned,mpfr_set_ui);
        _GMPFR_CONSTRUCTOR_FROM_TYPE(unsigned long,mpfr_set_ui);
        _GMPFR_CONSTRUCTOR_FROM_TYPE(double,mpfr_set_d);
        _GMPFR_CONSTRUCTOR_FROM_TYPE(long double,mpfr_set_ld);
        _GMPFR_CONSTRUCTOR_FROM_OBJECT(Gmpz,mpz(),mpfr_set_z);
        _GMPFR_CONSTRUCTOR_FROM_OBJECT(Gmpq,mpq(),mpfr_set_q);

#undef _GMPFR_CONSTRUCTOR_FROM_TYPE
#undef _GMPFR_CONSTRUCTOR_FROM_OBJECT

        // When Gmpfr is refence counted, we inherit the assignment
        // operator and the copy constructor from Handle_for.
#ifdef CGAL_GMPFR_NO_REFCOUNT
        Gmpfr& operator=(const Gmpfr &a){
                if(get_prec()!=a.get_prec())
                        set_prec(a.get_prec());
                mpfr_set(fr(),a.fr(),get_defrnd());
                return *this;
        }

        Gmpfr(const Gmpfr &a){
                mpfr_init2(fr(),a.get_prec());
                mpfr_set(fr(),a.fr(),GMP_RNDN);
        }
#endif

        Gmpfr(const Gmpfr &a,mp_rnd_t r,mp_prec_t p=get_defprec()){
#ifndef CGAL_GMPFR_NO_REFCOUNT
                if(p==a.get_prec()){
                        Gmpfr temp(a);
                        dont_clear_on_destruction();
                        swap(temp);
                }else
#endif
                {
                        mpfr_init2(fr(),p);
                        mpfr_set(fr(),a.fr(),r);
                }
        }

        Gmpfr(const Gmpfr &a,mp_prec_t p){
#ifndef CGAL_GMPFR_NO_REFCOUNT
                if(p==a.get_prec()){
                        Gmpfr temp(a);
                        dont_clear_on_destruction();
                        swap(temp);
                }else
#endif
                {
                        mpfr_init2(fr(),p);
                        mpfr_set(fr(),a.fr(),get_defrnd());
                }
        }

        // default rounding mode

        static mp_rnd_t get_defrnd();
        static void set_defrnd(mp_rnd_t);

        // default precision

        static mp_prec_t get_defprec();
        static void set_defprec(mp_prec_t);

        // precision of a single Gmpfr object

        mp_prec_t get_prec()const;
        void set_prec(mp_prec_t);
        void prec_round(mp_prec_t,mp_rnd_t);

        // mpfr global inexact flags

        static void clear_flags();
        static bool underflow_flag();
        static bool overflow_flag();
        static bool nan_flag();
        static bool inex_flag();
        static bool erange_flag();

        // arithmetics

        Gmpfr operator+()const;
        Gmpfr operator-()const;
        Gmpfr& operator+=(const Gmpfr&);
        Gmpfr& operator-=(const Gmpfr&);
        Gmpfr& operator*=(const Gmpfr&);
        Gmpfr& operator/=(const Gmpfr&);
        Gmpfr& operator%=(const Gmpfr&);
        Gmpfr& operator+=(long);
        Gmpfr& operator-=(long);
        Gmpfr& operator*=(long);
        Gmpfr& operator/=(long);
        Gmpfr& operator+=(unsigned long);
        Gmpfr& operator-=(unsigned long);
        Gmpfr& operator*=(unsigned long);
        Gmpfr& operator/=(unsigned long);
        Gmpfr& operator+=(int);
        Gmpfr& operator-=(int);
        Gmpfr& operator*=(int);
        Gmpfr& operator/=(int);
        Gmpfr& operator+=(const Gmpz&);
        Gmpfr& operator-=(const Gmpz&);
        Gmpfr& operator*=(const Gmpz&);
        Gmpfr& operator/=(const Gmpz&);
        Gmpfr& operator+=(const Gmpq&);
        Gmpfr& operator-=(const Gmpq&);
        Gmpfr& operator*=(const Gmpq&);
        Gmpfr& operator/=(const Gmpq&);

        bool is_zero()const;
        bool is_one()const;
        bool is_nan()const;
        bool is_inf()const;
        bool is_number()const;
        Sign sign()const;
        Gmpfr abs(mp_prec_t=0)const;
        Gmpfr sqrt(mp_prec_t=0)const;
        Gmpfr kthroot(int k,mp_prec_t=0)const;
        Gmpfr square(mp_prec_t=0)const;
        Comparison_result compare(const Gmpfr&)const;
        double to_double()const;
        std::pair<double,double>to_interval()const;
        std::pair<double,long>to_double_exp()const;
        std::pair<std::pair<double,double>,long>to_interval_exp()const;
        std::pair<Gmpz,long>to_integer_exp()const;
};




// --------------
// implementation
// --------------

// default rounding mode, handled by mpfr
inline
mp_rnd_t Gmpfr::get_defrnd(){
        return mpfr_get_default_rounding_mode();
};

inline
void Gmpfr::set_defrnd(mp_rnd_t rnd_mode){
        mpfr_set_default_rounding_mode(rnd_mode);
};

// default precision, handled by mpfr
inline
mp_prec_t Gmpfr::get_defprec(){
        return mpfr_get_default_prec();
};

inline
void Gmpfr::set_defprec(mp_prec_t prec){
        CGAL_assertion(prec>=MPFR_PREC_MIN&&prec<=MPFR_PREC_MAX);
        mpfr_set_default_prec(prec);
};

// precision of a single Gmpfr object

inline
mp_prec_t Gmpfr::get_prec()const{
        return mpfr_get_prec(fr());
}

inline
void Gmpfr::set_prec(mp_prec_t p){
        CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
#ifdef CGAL_GMPFR_NO_REFCOUNT
        if(p!=get_prec())
                mpfr_set_prec(fr(),p);
#else
        if(p!=get_prec()){
                if(unique()){
                        mpfr_set_prec(fr(),p);
                }else{
                        Gmpfr result(0,p);
                        mpfr_set_nan(result.fr());
                        swap(result);
                }
        }
#endif
}

inline
void Gmpfr::prec_round(mp_prec_t p,mp_rnd_t r){
        CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
#ifdef CGAL_GMPFR_NO_REFCOUNT
        if(p!=get_prec())
                mpfr_prec_round(fr(),p,r);
#else
        if(p!=get_prec()){
                if(unique()){
                        mpfr_prec_round(fr(),p,r);
                }else{
                        Gmpfr result(0,p);
                        mpfr_set(result.fr(),fr(),r);
                        swap(result);
                }
        }
#endif
}

// mpfr global inexact flags

inline
void Gmpfr::clear_flags(){
        mpfr_clear_flags();
}

inline
bool Gmpfr::underflow_flag(){
        return mpfr_underflow_p();
}

inline
bool Gmpfr::overflow_flag(){
        return mpfr_overflow_p();
}

inline
bool Gmpfr::nan_flag(){
        return mpfr_nanflag_p();
}

inline
bool Gmpfr::inex_flag(){
        return mpfr_inexflag_p();
}

inline
bool Gmpfr::erange_flag(){
        return mpfr_erangeflag_p();
}

// arithmetics

inline
Gmpfr Gmpfr::operator+()const{
        return(*this);
}

inline
Gmpfr Gmpfr::operator-()const{
        Gmpfr result(0,get_prec());
        mpfr_neg(result.fr(),fr(),GMP_RNDN);
        return result;
}

// _GMPFR_MAX_PREC returns the precision to be used to operate between
// *this and a number of another type or class. Currently, the maximum of
// *this' precision and the default precision is returned.
#define _GMPFR_MAX_PREC() (get_prec()>get_defprec()?get_prec():get_defprec())

// _GMPFR_MAX_PREC_2 returns the precision for the operation between Gmpfr
// objects *this and _b. Currently, it is the maximum of the precisions of
// *this and _b and the default precision.
// TODO: maybe we can rewrite this define optimally, maybe with an inline
#define _GMPFR_MAX_PREC_2(_b) \
        ( get_prec() >= mpfr_get_prec(_b.fr()) ? \
                ( get_prec()>(get_defprec())? \
                        get_prec():(get_defprec())): \
                ( mpfr_get_prec(_b.fr())>(get_defprec())? \
                        mpfr_get_prec(_b.fr()):(get_defprec())) \
        )

// _GMPFR_OBJECT_BINARY_OPERATOR defines an overloaded binary operator of
// the Gmpfr class, where the second parameter of the operator is an
// object. It behaves differently when the Gmpfr class is reference-counted
// or not.
#ifdef CGAL_GMPFR_NO_REFCOUNT
#define _GMPFR_OBJECT_BINARY_OPERATOR(_op,_class,_member,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(const _class &b){ \
                if(get_prec()>=get_defprec()) { \
                        _fun(fr(),fr(),b._member,get_defrnd()); \
                }else{ \
                        Gmpfr _temp(0,get_defprec()); \
                        _fun(_temp.fr(),fr(),b._member,get_defrnd()); \
                        mpfr_swap(_temp.fr(),fr()); \
                } \
                return *this; \
        }
#else
#define _GMPFR_OBJECT_BINARY_OPERATOR(_op,_class,_member,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(const _class &b){ \
                if(unique()){ \
                        if(get_prec()>get_defprec()) { \
                                _fun(fr(),fr(),b._member,get_defrnd()); \
                        }else{ \
                                Gmpfr _temp(0,get_defprec()); \
                                _fun(_temp.fr(),fr(),b._member,get_defrnd()); \
                                swap(_temp); \
                        } \
                }else{ \
                        Gmpfr result(0,_GMPFR_MAX_PREC()); \
                        _fun(result.fr(),fr(),b._member,get_defrnd()); \
                        swap(result); \
                } \
                return *this; \
        }
#endif

// _GMPFR_GMPFR_BINARY_OPERATOR is analogous to
// _GMPFR_OBJECT_BINARY_OPERATOR, and it is used when the second operand is
// another Gmpfr. The difference is the computation of the operation
// precision.
#ifdef CGAL_GMPFR_NO_REFCOUNT
#define _GMPFR_GMPFR_BINARY_OPERATOR(_op,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(const Gmpfr &b){ \
                mp_prec_t _p=_GMPFR_MAX_PREC_2(b); \
                if(_p==get_prec()) { \
                        _fun(fr(),fr(),b.fr(),get_defrnd()); \
                }else{ \
                        Gmpfr _temp(0,_p); \
                        _fun(_temp.fr(),fr(),b.fr(),get_defrnd()); \
                        mpfr_swap(_temp.fr(),fr()); \
                } \
                return *this; \
        }
#else
#define _GMPFR_GMPFR_BINARY_OPERATOR(_op,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(const Gmpfr &b){ \
                mp_prec_t _p=_GMPFR_MAX_PREC_2(b); \
                if(unique()&&(_p==get_prec())){ \
                        _fun(fr(),fr(),b.fr(),get_defrnd()); \
                }else{ \
                        Gmpfr result(0,_p); \
                        _fun(result.fr(),fr(),b.fr(),get_defrnd()); \
                        swap(result); \
                } \
                return *this; \
        }
#endif

// _GMPFR_TYPE_BINARY_OPERATOR is analogous to the
// _GMPFR_OBJECT_BINARY_OPERATOR, where the second parameter is a type
// instead of an object.
#ifdef CGAL_GMPFR_NO_REFCOUNT
#define _GMPFR_TYPE_BINARY_OPERATOR(_op,_type,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(_type x){ \
                if(get_prec()>=get_defprec()) { \
                        _fun(fr(),fr(),x,get_defrnd()); \
                }else{ \
                        Gmpfr _temp(0,get_defprec()); \
                        _fun(_temp.fr(),fr(),x,get_defrnd()); \
                        mpfr_swap(_temp.fr(),fr()); \
                } \
                return *this; \
        }
#else
#define _GMPFR_TYPE_BINARY_OPERATOR(_op,_type,_fun) \
        inline \
        Gmpfr& Gmpfr::_op(_type x){ \
                if(unique()){ \
                        if(get_prec()>get_defprec()) { \
                                _fun(fr(),fr(),x,get_defrnd()); \
                        }else{ \
                                Gmpfr _temp(0,get_defprec()); \
                                _fun(_temp.fr(),fr(),x,get_defrnd()); \
                                swap(_temp); \
                        } \
                }else{ \
                        Gmpfr result(0,_GMPFR_MAX_PREC()); \
                        _fun(result.fr(),fr(),x,get_defrnd()); \
                        swap(result); \
                } \
                return *this; \
        }
#endif

_GMPFR_GMPFR_BINARY_OPERATOR(operator+=,mpfr_add)
_GMPFR_GMPFR_BINARY_OPERATOR(operator-=,mpfr_sub)
_GMPFR_GMPFR_BINARY_OPERATOR(operator*=,mpfr_mul)
_GMPFR_GMPFR_BINARY_OPERATOR(operator/=,mpfr_div)
_GMPFR_GMPFR_BINARY_OPERATOR(operator%=,mpfr_div)

_GMPFR_TYPE_BINARY_OPERATOR(operator+=,long,mpfr_add_si)
_GMPFR_TYPE_BINARY_OPERATOR(operator-=,long,mpfr_sub_si)
_GMPFR_TYPE_BINARY_OPERATOR(operator*=,long,mpfr_mul_si)
_GMPFR_TYPE_BINARY_OPERATOR(operator/=,long,mpfr_div_si)

_GMPFR_TYPE_BINARY_OPERATOR(operator+=,unsigned long,mpfr_add_ui)
_GMPFR_TYPE_BINARY_OPERATOR(operator-=,unsigned long,mpfr_sub_ui)
_GMPFR_TYPE_BINARY_OPERATOR(operator*=,unsigned long,mpfr_mul_ui)
_GMPFR_TYPE_BINARY_OPERATOR(operator/=,unsigned long,mpfr_div_ui)

_GMPFR_TYPE_BINARY_OPERATOR(operator+=,int,mpfr_add_si)
_GMPFR_TYPE_BINARY_OPERATOR(operator-=,int,mpfr_sub_si)
_GMPFR_TYPE_BINARY_OPERATOR(operator*=,int,mpfr_mul_si)
_GMPFR_TYPE_BINARY_OPERATOR(operator/=,int,mpfr_div_si)

_GMPFR_OBJECT_BINARY_OPERATOR(operator+=,Gmpz,mpz(),mpfr_add_z)
_GMPFR_OBJECT_BINARY_OPERATOR(operator-=,Gmpz,mpz(),mpfr_sub_z)
_GMPFR_OBJECT_BINARY_OPERATOR(operator*=,Gmpz,mpz(),mpfr_mul_z)
_GMPFR_OBJECT_BINARY_OPERATOR(operator/=,Gmpz,mpz(),mpfr_div_z)

_GMPFR_OBJECT_BINARY_OPERATOR(operator+=,Gmpq,mpq(),mpfr_add_q)
_GMPFR_OBJECT_BINARY_OPERATOR(operator-=,Gmpq,mpq(),mpfr_sub_q)
_GMPFR_OBJECT_BINARY_OPERATOR(operator*=,Gmpq,mpq(),mpfr_mul_q)
_GMPFR_OBJECT_BINARY_OPERATOR(operator/=,Gmpq,mpq(),mpfr_div_q)

#undef _GMPFR_OBJECT_BINARY_OPERATOR
#undef _GMPFR_GMPFR_BINARY_OPERATOR
#undef _GMPFR_TYPE_BINARY_OPERATOR
#undef _GMPFR_MAX_PREC
#undef _GMPFR_MAX_PREC_2

inline
bool Gmpfr::is_zero()const{
        return mpfr_zero_p(fr());
}

inline
bool Gmpfr::is_one()const{
        return !mpfr_cmp_ui(fr(),1);
}

inline
bool Gmpfr::is_nan()const{
        return mpfr_nan_p(fr());
}

inline
bool Gmpfr::is_inf()const{
        return mpfr_inf_p(fr());
}

inline
bool Gmpfr::is_number()const{
        return mpfr_number_p(fr());
}

inline
Sign Gmpfr::sign()const{
        int s=mpfr_sgn(fr());
        return(s?(s>0?POSITIVE:NEGATIVE):ZERO);
}

inline
Gmpfr Gmpfr::abs(mp_prec_t p)const{
        if(!p)
                p=get_defprec();
        Gmpfr result(0,p);
        mpfr_abs(result.fr(),fr(),get_defrnd());
        return result;
}

inline
Gmpfr Gmpfr::sqrt(mp_prec_t p)const{
        if(!p)
                p=get_defprec();
        Gmpfr result(0,p);
        mpfr_sqrt(result.fr(),fr(),get_defrnd());
        return result;
}

inline
Gmpfr Gmpfr::kthroot(int k,mp_prec_t p)const{
        if(!p)
                p=get_defprec();
        Gmpfr result(0,p);
        if(k==3)
                mpfr_cbrt(result.fr(),fr(),get_defrnd());
        else
                mpfr_root(result.fr(),fr(),k,get_defrnd());
        return result;
}

inline
Gmpfr Gmpfr::square(mp_prec_t p)const{
        if(!p)
                p=get_defprec();
        Gmpfr result(0,p);
        mpfr_sqr(result.fr(),fr(),get_defrnd());
        return result;
}

inline
Comparison_result Gmpfr::compare(const Gmpfr& b)const{
        int c=mpfr_cmp(fr(),b.fr());
        return(c?(c>0?LARGER:SMALLER):EQUAL);
}

inline
double Gmpfr::to_double()const{
        return mpfr_get_d(fr(),get_defrnd());
}

// internal functions

inline
std::pair<double,long> Gmpfr::to_double_exp()const{
        long *e;
        double d=mpfr_get_d_2exp(e,fr(),get_defrnd());
        return std::make_pair(d,*e);
}

inline
std::pair<double,double>Gmpfr::to_interval()const{
        return std::make_pair(
                        mpfr_get_d(fr(),GMP_RNDD),
                        mpfr_get_d(fr(),GMP_RNDU));
}

inline
std::pair<std::pair<double,double>,long> Gmpfr::to_interval_exp()const{
        long *e1,*e2;
        double d_low=mpfr_get_d_2exp(e1,fr(),GMP_RNDD);
        double d_upp=mpfr_get_d_2exp(e2,fr(),GMP_RNDU);
        CGAL_assertion(*e1==*e2);
        return std::make_pair(std::make_pair(d_low,d_upp),*e1);
}

inline
std::pair<Gmpz,long> Gmpfr::to_integer_exp()const{
        Gmpz z;
        long e=mpfr_get_z_exp(z.mpz(),fr());
        return std::make_pair(z,e);
}

// input/output

inline
std::istream& operator>>(std::istream& is,Gmpfr &f){
        std::string s;
        is>>s;
        mpfr_set_str(f.fr(),s.c_str(),10,Gmpfr::get_defrnd());
        return is;
}

inline
std::ostream& operator<<(std::ostream& os,const Gmpfr &a){
        char *str;
        std::string s;
        mp_exp_t expptr;
        if(a.is_nan())
                return os<<"nan";
        if(a.is_inf())
                return os<<(a<0?"-inf":"+inf");
        str=mpfr_get_str(NULL,&expptr,10,0,a.fr(),Gmpfr::get_defrnd());
        s=(const char*)str;
        if(str[0]=='-'){
                os<<"-."<<s.substr(1);
        }else{
                os<<"."<<s;
        }
        os<<"e"<<expptr;
        mpfr_free_str(str);
        return os;
}

inline
std::ostream& print(std::ostream& os,const Gmpfr &a){
        mpz_t z;
        mpz_init(z);
        mp_exp_t e=mpfr_get_z_exp(z,a.fr());
        os<<Gmpz(z)<<"*2^"<<e;
        mpz_clear(z);
        return os;
}

// comparisons

inline
bool operator<(const Gmpfr &a,const Gmpfr &b){
        return mpfr_less_p(a.fr(),b.fr());
}

inline
bool operator==(const Gmpfr &a,const Gmpfr &b){
        return mpfr_equal_p(a.fr(),b.fr());
}

inline
bool operator<(const Gmpfr &a,long b){
        return(mpfr_cmp_si(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,long b){
        return(mpfr_cmp_si(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,long b){
        return !mpfr_cmp_si(a.fr(),b);
}

inline
bool operator<(const Gmpfr &a,unsigned long b){
        return(mpfr_cmp_ui(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,unsigned long b){
        return(mpfr_cmp_ui(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,unsigned long b){
        return !mpfr_cmp_ui(a.fr(),b);
}

inline
bool operator<(const Gmpfr &a,int b){
        return(mpfr_cmp_si(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,int b){
        return(mpfr_cmp_si(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,int b){
        return !mpfr_cmp_si(a.fr(),b);
}

inline
bool operator<(const Gmpfr &a,double b){
        return(mpfr_cmp_d(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,double b){
        return(mpfr_cmp_d(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,double b){
        return !mpfr_cmp_d(a.fr(),b);
}

inline
bool operator<(const Gmpfr &a,long double b){
        return(mpfr_cmp_ld(a.fr(),b)<0);
}

inline
bool operator>(const Gmpfr &a,long double b){
        return(mpfr_cmp_ld(a.fr(),b)>0);
}

inline
bool operator==(const Gmpfr &a,long double b){
        return !mpfr_cmp_ld(a.fr(),b);
}

inline
bool operator<(const Gmpfr &a,const Gmpz &b){
        return(mpfr_cmp_z(a.fr(),b.mpz())<0);
}

inline
bool operator>(const Gmpfr &a,const Gmpz &b){
        return(mpfr_cmp_z(a.fr(),b.mpz())>0);
}

inline
bool operator==(const Gmpfr &a,const Gmpz &b){
        return !mpfr_cmp_z(a.fr(),b.mpz());
}

inline
bool operator<(const Gmpfr &a,const Gmpq &b){
        return(mpfr_cmp_q(a.fr(),b.mpq())<0);
}

inline
bool operator>(const Gmpfr &a,const Gmpq &b){
        return(mpfr_cmp_q(a.fr(),b.mpq())>0);
}

inline
bool operator==(const Gmpfr &a,const Gmpq &b){
        return !mpfr_cmp_q(a.fr(),b.mpq());
}

inline
Gmpfr min BOOST_PREVENT_MACRO_SUBSTITUTION(const Gmpfr& x,const Gmpfr& y){
        return (x<=y)?x:y;
}

inline
Gmpfr max BOOST_PREVENT_MACRO_SUBSTITUTION(const Gmpfr& x,const Gmpfr& y){
        return (x>=y)?x:y;
}

} // namespace CGAL 

#endif  // CGAL_GMPFR_TYPE_H

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab

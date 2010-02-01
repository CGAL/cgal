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

#ifndef CGAL_GMPFI_TYPE_H
#define CGAL_GMPFI_TYPE_H

#include <CGAL/basic.h>
#include <gmp.h>
#include <mpfr.h>
#include <CGAL/GMP/Gmpfr_type.h>
#include <mpfi.h>
#include <boost/operators.hpp>
#include <CGAL/Handle_for.h>
#include <CGAL/Uncertain.h>
#ifdef CGAL_HAS_THREADS
#  include <boost/thread/tss.hpp>
#endif
#include <limits>

namespace CGAL{

class Gmpfi;

Uncertain<bool> operator<(const Gmpfi&,const Gmpfi&);
Uncertain<bool> operator==(const Gmpfi&,const Gmpfi&);

Uncertain<bool> operator<(const Gmpfi&,const Gmpfr&);
Uncertain<bool> operator>(const Gmpfi&,const Gmpfr&);
Uncertain<bool> operator==(const Gmpfi&,const Gmpfr&);

Uncertain<bool> operator<(const Gmpfi&,long);
Uncertain<bool> operator>(const Gmpfi&,long);
Uncertain<bool> operator==(const Gmpfi&,long);

Uncertain<bool> operator<(const Gmpfi&,unsigned long);
Uncertain<bool> operator>(const Gmpfi&,unsigned long);
Uncertain<bool> operator==(const Gmpfi&,unsigned long);

Uncertain<bool> operator<(const Gmpfi&,int);
Uncertain<bool> operator>(const Gmpfi&,int);
Uncertain<bool> operator==(const Gmpfi&,int);

Uncertain<bool> operator<(const Gmpfi&,double);
Uncertain<bool> operator>(const Gmpfi&,double);
Uncertain<bool> operator==(const Gmpfi&,double);

Uncertain<bool> operator<(const Gmpfi&,long double);
Uncertain<bool> operator>(const Gmpfi&,long double);
Uncertain<bool> operator==(const Gmpfi&,long double);

Uncertain<bool> operator<(const Gmpfi&,const Gmpz&);
Uncertain<bool> operator>(const Gmpfi&,const Gmpz&);
Uncertain<bool> operator==(const Gmpfi&,const Gmpz&);

Uncertain<bool> operator<(const Gmpfi&,const Gmpq&);
Uncertain<bool> operator>(const Gmpfi&,const Gmpq&);
Uncertain<bool> operator==(const Gmpfi&,const Gmpq&);

struct Gmpfi_rep{
        mpfi_t floating_point_interval;
        bool clear_on_destruction;
        Gmpfi_rep():clear_on_destruction(true){}
        ~Gmpfi_rep(){
                if(clear_on_destruction)
                        mpfi_clear(floating_point_interval);
        }
};

// the default precision of Gmpfi is the size of a double's mantissa
#ifdef IEEE_DBL_MANT_DIG
#  define _GMPFI_DEFAULT_PRECISION IEEE_DBL_MANT_DIG
#else
#  define _GMPFI_DEFAULT_PRECISION 53
#endif

// the default precision is a variable local to each thread in multithreaded
// environments, or a global variable otherwise
#ifdef CGAL_HAS_THREADS
        boost::thread_specific_ptr<mp_prec_t> Gmpfi_default_precision_;
#else
        mp_prec_t Gmpfi_default_precision=_GMPFI_DEFAULT_PRECISION;
#endif

class Gmpfi:
        Handle_for<Gmpfi_rep>,
        boost::ordered_euclidian_ring_operators1<Gmpfi,
        boost::ordered_euclidian_ring_operators2<Gmpfi,Gmpfr,
        boost::ordered_euclidian_ring_operators2<Gmpfi,long,
        boost::ordered_euclidian_ring_operators2<Gmpfi,unsigned long,
        boost::ordered_euclidian_ring_operators2<Gmpfi,int,
        boost::totally_ordered2<Gmpfi,double,
        boost::totally_ordered2<Gmpfi,long double,
        boost::ordered_euclidian_ring_operators2<Gmpfi,Gmpz,
        boost::ordered_euclidian_ring_operators2<Gmpfi,Gmpq
        > > > > > > > > >
{
        typedef Handle_for<Gmpfi_rep>   Base;

        public:

        typedef Gmpfr::Precision_type   Precision_type;

        // access

        inline mpfi_srcptr mpfi()const{
                return Ptr()->floating_point_interval;
        }

        inline mpfr_srcptr left_mpfr()const{
                return &(mpfi()->left);
        }

        inline mpfr_srcptr right_mpfr()const{
                return &(mpfi()->right);
        }

        inline mpfi_ptr mpfi(){
                return ptr()->floating_point_interval;
        }

        inline Gmpfr inf()const{
                return Gmpfr(left_mpfr());
        }

        inline Gmpfr sup()const{
                return Gmpfr(right_mpfr());
        }

        inline
        void dont_clear_on_destruction(){
                ptr()->clear_on_destruction=false;
        }

        // construction

        Gmpfi(){
                mpfi_init(mpfi());
        }

#define _GMPFI_CONSTRUCTOR_FROM_SCALAR(_type) \
        Gmpfi(const _type &t, \
              Gmpfi::Precision_type p=Gmpfi::get_default_precision()){ \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfr l(t,std::round_toward_neg_infinity,p), \
                      r(t,std::round_toward_infinity,p); \
                l.dont_clear_on_destruction(); \
                r.dont_clear_on_destruction(); \
                mpfi()->left=*(l.fr()); \
                mpfi()->right=*(r.fr()); \
        }

_GMPFI_CONSTRUCTOR_FROM_SCALAR(long);
_GMPFI_CONSTRUCTOR_FROM_SCALAR(unsigned long);
_GMPFI_CONSTRUCTOR_FROM_SCALAR(int);
_GMPFI_CONSTRUCTOR_FROM_SCALAR(double);
_GMPFI_CONSTRUCTOR_FROM_SCALAR(long double);
_GMPFI_CONSTRUCTOR_FROM_SCALAR(Gmpz);
_GMPFI_CONSTRUCTOR_FROM_SCALAR(Gmpq);

        Gmpfi(mpfi_srcptr i,Gmpfi::Precision_type p=0){
                if((p==0)||
                   (p==mpfr_get_prec(&(i->left))
                    &&p==mpfr_get_prec(&(i->right)))){
                        mpfi()->left=i->left;
                        mpfi()->right=i->right;
                        dont_clear_on_destruction();
                }else{
                        CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                        mpfi_init2(mpfi(),p);
                        mpfi_set(mpfi(),i);
                }
        }

        Gmpfi(const Gmpfr &f,
              Gmpfi::Precision_type p=Gmpfi::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                mpfi_init2(mpfi(),p);
                mpfi_set_fr(mpfi(),f.fr());
        }

        Gmpfi(std::pair<const Gmpfr,const Gmpfr> endpoints,
              Gmpfi::Precision_type p=Gmpfi::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                mpfi_init2(mpfi(),p);
                mpfi_interv_fr(
                                mpfi(),
                                endpoints.first.fr(),
                                endpoints.second.fr());
        }

        template<class L,class R>
        Gmpfi(std::pair<const L&,const R&> endpoints,
              Gmpfi::Precision_type p=get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                Gmpfr l(endpoints.first,std::round_toward_neg_infinity,p),
                      r(endpoints.second,std::round_toward_infinity,p);
                l.dont_clear_on_destruction();
                r.dont_clear_on_destruction();
                mpfi()->left=*(l.fr());
                mpfi()->right=*(r.fr());
        }

        // copy constructor
        Gmpfi(const Gmpfi &a,Gmpfi::Precision_type p){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                if(p==mpfr_get_prec(a.left_mpfr())&&
                                p==mpfr_get_prec(a.right_mpfr())){
                        Gmpfi temp(a);
                        dont_clear_on_destruction();
                        swap(temp);
                }else{
                        mpfi_init2(mpfi(),p);
                        mpfi_set(mpfi(),a.mpfi());
                }
        }

        // default precision

#ifdef CGAL_HAS_THREADS
        static void init_precision_for_thread();
#endif
        static inline Gmpfi::Precision_type get_default_precision();
        static inline Gmpfi::Precision_type set_default_precision(
                                                Gmpfi::Precision_type prec);

        // precision of a single Gmpfi object

        Gmpfi::Precision_type get_precision()const;
        Gmpfi round(Gmpfi::Precision_type)const;

        // arithmetics

        Gmpfi operator+()const;
        Gmpfi operator-()const;

#define _GMPFI_DECLARE_OPERATORS(_type) \
        Gmpfi& operator+=(_type); \
        Gmpfi& operator-=(_type); \
        Gmpfi& operator*=(_type); \
        Gmpfi& operator/=(_type);

        _GMPFI_DECLARE_OPERATORS(const Gmpfi&)
        _GMPFI_DECLARE_OPERATORS(const Gmpfr&)
        _GMPFI_DECLARE_OPERATORS(long)
        _GMPFI_DECLARE_OPERATORS(unsigned long)
        _GMPFI_DECLARE_OPERATORS(int)
        _GMPFI_DECLARE_OPERATORS(const Gmpz&)
        _GMPFI_DECLARE_OPERATORS(const Gmpq&)

#undef _GMPFI_DECLARE_OPERATORS

#define _GMPFI_DECLARE_STATIC_FUNCTION(_f,_t1,_t2) \
        static Gmpfi _f (_t1,_t2,Gmpfi::Precision_type=0);

#define _GMPFI_DECLARE_STATIC_FUNCTIONS(_t1,_t2) \
        _GMPFI_DECLARE_STATIC_FUNCTION(add,_t1,_t2) \
        _GMPFI_DECLARE_STATIC_FUNCTION(sub,_t1,_t2) \
        _GMPFI_DECLARE_STATIC_FUNCTION(mul,_t1,_t2) \
        _GMPFI_DECLARE_STATIC_FUNCTION(div,_t1,_t2)

#define _GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(_t) \
        _GMPFI_DECLARE_STATIC_FUNCTIONS(const Gmpfi&,_t) \
        _GMPFI_DECLARE_STATIC_FUNCTIONS(_t,const Gmpfi&)

        _GMPFI_DECLARE_STATIC_FUNCTIONS(const Gmpfi&,const Gmpfi&)
        _GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(const Gmpfr&)
        _GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(long)
        _GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(unsigned long)
        _GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(int)
        _GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(const Gmpz&)
        _GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(const Gmpq&)

#undef _GMPFI_DECLARE_STATIC_FUNCTION

        Gmpfi abs(Gmpfi::Precision_type=Gmpfi::get_default_precision())const;
        Gmpfi sqrt(Gmpfi::Precision_type=Gmpfi::get_default_precision())const;
        Gmpfi cbrt(Gmpfi::Precision_type=Gmpfi::get_default_precision())const;
        Gmpfi kthroot(int,
                      Gmpfi::Precision_type=Gmpfi::get_default_precision()
                     )const;
        Gmpfi square(Gmpfi::Precision_type=Gmpfi::get_default_precision())const;

        // comparison and query functions

        bool is_point()const;
        bool is_same(const Gmpfi&)const;
        bool do_overlap(const Gmpfi&)const;
        Uncertain<bool> is_zero()const;
        Uncertain<bool> is_one()const;
        bool is_nan()const;
        bool is_inf()const;
        bool is_number()const;
        Uncertain<Sign> sign()const;
        Uncertain<bool> is_positive()const;
        Uncertain<bool> is_negative()const;
        Uncertain<bool> is_square()const;
        Uncertain<bool> is_square(Gmpfi&)const;
        Uncertain<bool> divides(const Gmpfi&,
                                Gmpfi&,
                                Gmpfi::Precision_type=
                                        Gmpfi::get_default_precision()
                               )const;
        Uncertain<Comparison_result> compare(const Gmpfi&)const;

        // conversion functions

        double to_double()const;
        std::pair<double,double> to_interval()const;
        std::pair<double,long> to_double_exp()const;
        std::pair<std::pair<double,double>,long> to_interval_exp()const;
};




// --------------
// implementation
// --------------

// default precision
#ifdef CGAL_HAS_THREADS
void Gmpfi::init_precision_for_thread(){
        CGAL_precondition(Gmpfi_default_precision_.get()==NULL);
        Gmpfi_default_precision_.reset(
                new mp_prec_t(_GMPFI_DEFAULT_PRECISION));
}
#endif

inline
Gmpfi::Precision_type Gmpfi::get_default_precision(){
#ifdef CGAL_HAS_THREADS
        if(Gmpfi_default_precision_.get()==NULL)
                Gmpfi::init_precision_for_thread();
        return *Gmpfi_default_precision_.get();
#else
        return Gmpfi_default_precision;
#endif
}

inline
Gmpfi::Precision_type Gmpfi::set_default_precision(Gmpfi::Precision_type prec){
        Gmpfi::Precision_type old_prec=Gmpfi::get_default_precision();
        CGAL_assertion(prec>=MPFR_PREC_MIN&&prec<=MPFR_PREC_MAX);
#ifdef CGAL_HAS_THREADS
        *Gmpfi_default_precision_.get()=prec;
#else
        Gmpfi_default_precision=prec;
#endif
        return old_prec;
}

// precision of a single Gmpfi object

inline
Gmpfi::Precision_type Gmpfi::get_precision()const{
        return mpfi_get_prec(mpfi());
}

inline
Gmpfi Gmpfi::round(Gmpfi::Precision_type p)const{
        CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
        return Gmpfi(*this,p);
}

// arithmetics

inline
Gmpfi Gmpfi::operator+()const{
        return(*this);
}

inline
Gmpfi Gmpfi::operator-()const{
        Gmpfi result(0,get_precision());
        mpfi_neg(result.mpfi(),mpfi());
        return result;
}

// _GMPFI_BALANCE_ENDPOINTS checks if both endpoints of the interval have
// the same precision. If not, it rounds the one with the smallest
// precision.
#define _GMPFI_BALANCE_ENDPOINTS \
        if(mpfr_get_prec(left_mpfr())<mpfr_get_prec(right_mpfr())){\
                mpfr_round_prec(&(mpfi()->left), \
                                GMP_RNDD, \
                                mpfr_get_prec(right_mpfr())); \
        }else{ \
                if(mpfr_get_prec(left_mpfr())>mpfr_get_prec(right_mpfr())){\
                        mpfr_round_prec(&(mpfi()->right), \
                                        GMP_RNDU, \
                                        mpfr_get_prec(left_mpfr())); \
                } \
        }

// _GMPFI_OBJECT_BINARY_OPERATOR defines an overloaded binary operator of
// the Gmpfi class, where the operated object belongs to another class,
// which represents a point (as opposition to an interval). The operation
// will be performed using the biggest precision of the endpoints of this
// Gmpfi object. That means that if endpoints have different precision, one
// of them (the one with the biggest precision) will be rounded. This is
// not a problem when the object is not unique, since a new Gmpfi object
// will be created with the endpoints having the correct precision.
#define _GMPFI_OBJECT_BINARY_OPERATOR(_op,_class,_member,_fun) \
        inline \
        Gmpfi& Gmpfi::_op(const _class &b){ \
                if(unique()){ \
                        _GMPFI_BALANCE_ENDPOINTS \
                        _fun(mpfi(),mpfi(),b._member); \
                }else{ \
                        Gmpfi result(0,get_precision()); \
                        _fun(result.mpfi(),mpfi(),b._member); \
                        swap(result); \
                } \
                return(*this); \
        }

// _GMPFI_GMPFI_BINARY_OPERATOR defines an overloaded binary operator of
// the Gmpfi class, where the operated object is also a Gmpfi object.
// The operation will be performed using the biggest precision of the
// endpoints of both intervals. The endpoints of target object will be
// rounded accordingly before the operation.
#define _GMPFI_GMPFI_BINARY_OPERATOR(_op,_fun) \
        inline \
        Gmpfi& Gmpfi::_op(const Gmpfi &fi){ \
                if(unique()){ \
                        if(get_precision()<fi.get_precision()){ \
                                mpfi_round_prec(mpfi(),fi.get_precision()); \
                        }else{ \
                                _GMPFI_BALANCE_ENDPOINTS \
                        } \
                        _fun(mpfi(),mpfi(),fi.mpfi()); \
                }else{ \
                        Gmpfi result(0, \
                                     get_precision()<fi.get_precision()? \
                                        fi.get_precision(): \
                                        get_precision()); \
                        _fun(result.mpfi(),mpfi(),fi.mpfi()); \
                        swap(result); \
                } \
                return(*this); \
        }

// _GMPFI_TYPE_BINARY_OPERATOR defines an overloaded binary operator of
// the Gmpfi class, where the operated belongs to a c++ type. Precision of
// the operation is defined in the same manner that in
// _GMPFI_OBJECT_BINARY_OPERATOR.
#define _GMPFI_TYPE_BINARY_OPERATOR(_op,_type,_fun) \
        inline \
        Gmpfi& Gmpfi::_op(_type x){ \
                if(unique()){ \
                        _GMPFI_BALANCE_ENDPOINTS \
                        _fun(mpfi(),mpfi(),x); \
                }else{ \
                        Gmpfi result(0,get_precision()); \
                        _fun(result.mpfi(),mpfi(),x); \
                        swap(result); \
                } \
                return *this; \
        }

_GMPFI_GMPFI_BINARY_OPERATOR(operator+=,mpfi_add)
_GMPFI_GMPFI_BINARY_OPERATOR(operator-=,mpfi_sub)
_GMPFI_GMPFI_BINARY_OPERATOR(operator*=,mpfi_mul)
_GMPFI_GMPFI_BINARY_OPERATOR(operator/=,mpfi_div)

_GMPFI_OBJECT_BINARY_OPERATOR(operator+=,Gmpfr,fr(),mpfi_add_fr)
_GMPFI_OBJECT_BINARY_OPERATOR(operator-=,Gmpfr,fr(),mpfi_sub_fr)
_GMPFI_OBJECT_BINARY_OPERATOR(operator*=,Gmpfr,fr(),mpfi_mul_fr)
_GMPFI_OBJECT_BINARY_OPERATOR(operator/=,Gmpfr,fr(),mpfi_div_fr)

_GMPFI_TYPE_BINARY_OPERATOR(operator+=,long,mpfi_add_si)
_GMPFI_TYPE_BINARY_OPERATOR(operator-=,long,mpfi_sub_si)
_GMPFI_TYPE_BINARY_OPERATOR(operator*=,long,mpfi_mul_si)
_GMPFI_TYPE_BINARY_OPERATOR(operator/=,long,mpfi_div_si)

_GMPFI_TYPE_BINARY_OPERATOR(operator+=,unsigned long,mpfi_add_ui)
_GMPFI_TYPE_BINARY_OPERATOR(operator-=,unsigned long,mpfi_sub_ui)
_GMPFI_TYPE_BINARY_OPERATOR(operator*=,unsigned long,mpfi_mul_ui)
_GMPFI_TYPE_BINARY_OPERATOR(operator/=,unsigned long,mpfi_div_ui)

_GMPFI_TYPE_BINARY_OPERATOR(operator+=,int,mpfi_add_si)
_GMPFI_TYPE_BINARY_OPERATOR(operator-=,int,mpfi_sub_si)
_GMPFI_TYPE_BINARY_OPERATOR(operator*=,int,mpfi_mul_si)
_GMPFI_TYPE_BINARY_OPERATOR(operator/=,int,mpfi_div_si)

_GMPFI_OBJECT_BINARY_OPERATOR(operator+=,Gmpz,mpz(),mpfi_add_z)
_GMPFI_OBJECT_BINARY_OPERATOR(operator-=,Gmpz,mpz(),mpfi_sub_z)
_GMPFI_OBJECT_BINARY_OPERATOR(operator*=,Gmpz,mpz(),mpfi_mul_z)
_GMPFI_OBJECT_BINARY_OPERATOR(operator/=,Gmpz,mpz(),mpfi_div_z)

_GMPFI_OBJECT_BINARY_OPERATOR(operator+=,Gmpq,mpq(),mpfi_add_q)
_GMPFI_OBJECT_BINARY_OPERATOR(operator-=,Gmpq,mpq(),mpfi_sub_q)
_GMPFI_OBJECT_BINARY_OPERATOR(operator*=,Gmpq,mpq(),mpfi_mul_q)
_GMPFI_OBJECT_BINARY_OPERATOR(operator/=,Gmpq,mpq(),mpfi_div_q)

#undef _GMPFI_GMPFI_BINARY_OPERATOR
#undef _GMPFI_OBJECT_BINARY_OPERATOR
#undef _GMPFI_TYPE_BINARY_OPERATOR

// the static arithmetic functions are defined in a separate file
#include <CGAL/GMP/Gmpfi_type_static.h>

#define _GMPFI_ARITHMETIC_FUNCTION(_name,_fun) \
        inline \
        Gmpfi Gmpfi::_name (Gmpfi::Precision_type p)const{ \
                Gmpfi result(0,p); \
                _fun(result.mpfi(),mpfi()); \
                return result; \
        }

_GMPFI_ARITHMETIC_FUNCTION(abs,mpfi_abs)
_GMPFI_ARITHMETIC_FUNCTION(sqrt,mpfi_sqrt)

inline
Gmpfi Gmpfi::cbrt(Gmpfi::Precision_type p)const{
        // MPFI does not provide a cubic root function
        Gmpfi result(0,p);
        mpfr_cbrt(&(result.mpfi())->left,left_mpfr(),GMP_RNDD);
        mpfr_cbrt(&(result.mpfi())->right,right_mpfr(),GMP_RNDU);
        return result;
}

inline
Gmpfi Gmpfi::kthroot(int k,Gmpfi::Precision_type p)const{
        // MPFI does not provide k-th root functions
        Gmpfi result(0,p);
        mpfr_root(&(result.mpfi())->left,left_mpfr(),k,GMP_RNDD);
        mpfr_root(&(result.mpfi())->right,right_mpfr(),k,GMP_RNDU);
        return result;
}

_GMPFI_ARITHMETIC_FUNCTION(square,mpfi_sqr)

// comparison and query functions

inline
bool Gmpfi::is_point()const{
        return mpfr_equal_p(left_mpfr(),right_mpfr())!=0;
}

inline
bool Gmpfi::is_same(const Gmpfi &b)const{
        return(mpfr_equal_p(left_mpfr(),b.left_mpfr())!=0 &&
                mpfr_equal_p(right_mpfr(),b.right_mpfr())!=0);
}

inline
bool Gmpfi::do_overlap(const Gmpfi &b)const{
        if(mpfr_lessequal_p(left_mpfr(),b.left_mpfr())!=0)
                return mpfr_lessequal_p(b.left_mpfr(),right_mpfr())!=0;
        else
                return mpfr_lessequal_p(left_mpfr(),b.right_mpfr())!=0;
}

inline
Uncertain<bool> Gmpfi::is_zero()const{
        if(mpfr_zero_p(&mpfi()->left)!=0 && mpfr_zero_p(&mpfi()->right)!=0)
                return true;
        if(mpfi_has_zero(mpfi())!=0)
                return Uncertain<bool>::indeterminate();
        return false;
}

inline
Uncertain<bool> Gmpfi::is_one()const{
        if(mpfr_cmp_ui(left_mpfr(),1)==0 && mpfr_cmp_ui(right_mpfr(),1)==0)
                return true;
        if(mpfi_is_inside_ui(1,mpfi())!=0)
                return Uncertain<bool>::indeterminate();
        return false;
}

inline
bool Gmpfi::is_nan()const{
        return mpfi_nan_p(mpfi())!=0;
}

inline
bool Gmpfi::is_inf()const{
        return mpfi_inf_p(mpfi())!=0;
}

inline
bool Gmpfi::is_number()const{
        return mpfi_bounded_p(mpfi())!=0;
}

inline
Uncertain<Sign> Gmpfi::sign()const{
        int s_left=mpfr_sgn(left_mpfr());
        if(s_left>0)
                return POSITIVE;
        if(s_left==0){
                if(mpfr_zero_p(right_mpfr())!=0)
                        return ZERO;
                else
                        return Uncertain<Sign>::indeterminate();
        }
        if(mpfr_sgn(right_mpfr())<0)
                return NEGATIVE;
        return Uncertain<Sign>::indeterminate();
}

inline
Uncertain<bool> Gmpfi::is_positive()const{
        if(mpfr_sgn(left_mpfr())>0)
                return true;
        if(mpfr_sgn(right_mpfr())<=0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> Gmpfi::is_negative()const{
        if(mpfr_sgn(right_mpfr())<0)
                return true;
        if(mpfr_sgn(left_mpfr())>=0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> Gmpfi::is_square()const{
        if(mpfr_sgn(left_mpfr())>=0)
                return true;
        if(mpfr_sgn(right_mpfr())<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> Gmpfi::is_square(Gmpfi &y)const{
        if(mpfr_sgn(left_mpfr())>=0){
                y=sqrt();
                return true;
        }
        if(mpfr_sgn(right_mpfr())<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> Gmpfi::divides(const Gmpfi &n,Gmpfi &c,Gmpfi::Precision_type p
                              )const{
        if(mpfr_zero_p(&mpfi()->left)!=0 && mpfr_zero_p(&mpfi()->right)!=0)
                return false;
        if(mpfi_has_zero(mpfi())!=0)
                return Uncertain<bool>::indeterminate();
        c=Gmpfi::div(n,*this,p);
        return true;
}

inline
Uncertain<Comparison_result> Gmpfi::compare(const Gmpfi& b)const{
        if(mpfr_greater_p(left_mpfr(),b.right_mpfr())!=0)
                return LARGER;
        if(mpfr_greater_p(b.left_mpfr(),right_mpfr())!=0)
                return SMALLER;
        if(mpfr_equal_p(left_mpfr(),b.right_mpfr())!=0 &&
                        mpfr_equal_p(b.left_mpfr(),right_mpfr())!=0)
                return EQUAL;
        return Uncertain<Comparison_result>::indeterminate();
}

// conversion functions

inline
double Gmpfi::to_double()const{
        return mpfi_get_d(mpfi());
}

inline
std::pair<double,double> Gmpfi::to_interval()const{
        double d_low=mpfr_get_d(left_mpfr(),GMP_RNDD);
        double d_upp=mpfr_get_d(right_mpfr(),GMP_RNDU);
        CGAL_assertion(std::numeric_limits<double>::has_infinity);
        // if an endpoint is finite and its double is infinity, we overflow
        if(mpfr_inf_p(left_mpfr())==0&&
                        d_low==std::numeric_limits<double>::infinity())
                mpfr_set_underflow();
        if(mpfr_inf_p(right_mpfr())==0&&
                        d_upp==std::numeric_limits<double>::infinity())
                mpfr_set_overflow();
        return std::make_pair(d_low,d_upp);
}

inline
std::pair<double,long> Gmpfi::to_double_exp()const{
        mpfr_t middle;
        long *e;
        mpfr_init2(middle,53);
        mpfi_get_fr(middle,mpfi());
        double d=mpfr_get_d_2exp(e,middle,mpfr_get_default_rounding_mode());
        mpfr_clear(middle);
        return std::make_pair(d,*e);
}

inline
std::pair<std::pair<double,double>,long> Gmpfi::to_interval_exp()const{
        long *e1,*e2;
        double d_low=mpfr_get_d_2exp(e1,left_mpfr(),GMP_RNDD);
        double d_upp=mpfr_get_d_2exp(e2,right_mpfr(),GMP_RNDU);
        if(e1<e2)
                d_upp=d_upp/pow(2.,(double)((*e2)-(*e1)));
        else if(e1>e2){
                d_low=d_low/pow(2.,(double)((*e1)-(*e2)));
                *e1=*e2;
        }
        return std::make_pair(std::make_pair(d_low,d_upp),*e1);
}

// input/output

// This function reads an interval from the istream. It has the form
// [inf,sup], where each one of inf and sup is read as a Gmpfr. Then, they
// are rounded accordingly. The input may contain spaces between the
// brackets and the numbers and the numbers and the comma. The result is
// undefined when the input is malformed.
inline
std::istream& operator>>(std::istream& is,Gmpfi &f){
        Gmpfr left,right;
        std::istream::int_type c;
        std::ios::fmtflags old_flags = is.flags();
        is.unsetf(std::ios::skipws);
        gmpz_eat_white_space(is);
        c=is.get();
        if(c!='['){
                invalid_number:
                is.setstate(std::ios_base::failbit);
                is.flags(old_flags);
                return is;
        }
        gmpz_eat_white_space(is);
        is>>left;
        c=is.get();
        if(c!=',')
                goto invalid_number;
        is>>right;
        gmpz_eat_white_space(is);
        c=is.get();
        if(c!=']')
                goto invalid_number;
        Gmpfr::Precision_type p=left.get_precision()>right.get_precision()?
                                left.get_precision():
                                right.get_precision();
        f=Gmpfi(std::make_pair(left,right),(Gmpfi::Precision_type)p);
        return is;
}

inline
std::ostream& operator<<(std::ostream& os,const Gmpfi &a){
        return(os<<"["<<a.inf()<<","<<a.sup()<<"]");
}

// comparisons

inline
Uncertain<bool> operator<(const Gmpfi &a,const Gmpfi &b){
        int c=mpfi_cmp(a.mpfi(),b.mpfi());
        if(c<0)
                return true;
        if(c>0 || (a.is_point() && b.is_point()))
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator==(const Gmpfi &a,const Gmpfi &b){
        if(mpfr_less_p(a.right_mpfr(),b.left_mpfr()) ||
                        mpfr_less_p(b.right_mpfr(),a.left_mpfr()))
                return false;
        if(mpfr_equal_p(a.left_mpfr(),b.right_mpfr()) &&
                        mpfr_equal_p(b.left_mpfr(),a.right_mpfr()))
                return true;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator<(const Gmpfi &a,const Gmpfr &b){
        if(mpfr_cmp(a.right_mpfr(),b.fr())<0)
                return true;
        if(mpfr_cmp(a.left_mpfr(),b.fr())>0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator>(const Gmpfi &a,const Gmpfr &b){
        if(mpfr_cmp(a.left_mpfr(),b.fr())>0)
                return true;
        if(mpfr_cmp(a.right_mpfr(),b.fr())<0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator==(const Gmpfi &a,const Gmpfr &b){
        if(a.is_point())
                return(mpfr_cmp(a.left_mpfr(),b.fr())?false:true);
        if(mpfr_cmp(a.left_mpfr(),b.fr())>0 ||
                        mpfr_cmp(a.right_mpfr(),b.fr())<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator<(const Gmpfi &a,long b){
        if(mpfr_cmp_si(a.right_mpfr(),b)<0)
                return true;
        if(mpfr_cmp_si(a.left_mpfr(),b)>0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator>(const Gmpfi &a,long b){
        if(mpfr_cmp_si(a.left_mpfr(),b)>0)
                return true;
        if(mpfr_cmp_si(a.right_mpfr(),b)<0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator==(const Gmpfi &a,long b){
        if(a.is_point())
                return(mpfr_cmp_si(a.left_mpfr(),b)?false:true);
        if(mpfr_cmp_si(a.left_mpfr(),b)>0 || mpfr_cmp_si(a.right_mpfr(),b)<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator<(const Gmpfi &a,unsigned long b){
        if(mpfr_cmp_ui(a.right_mpfr(),b)<0)
                return true;
        if(mpfr_cmp_ui(a.left_mpfr(),b)>0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator>(const Gmpfi &a,unsigned long b){
        if(mpfr_cmp_ui(a.left_mpfr(),b)>0)
                return true;
        if(mpfr_cmp_ui(a.right_mpfr(),b)<0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator==(const Gmpfi &a,unsigned long b){
        if(a.is_point())
                return(mpfr_cmp_ui(a.left_mpfr(),b)?false:true);
        if(mpfr_cmp_ui(a.left_mpfr(),b)>0 || mpfr_cmp_ui(a.right_mpfr(),b)<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator<(const Gmpfi &a,int b){
        if(mpfr_cmp_si(a.right_mpfr(),b)<0)
                return true;
        if(mpfr_cmp_si(a.left_mpfr(),b)>0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator>(const Gmpfi &a,int b){
        if(mpfr_cmp_si(a.left_mpfr(),b)>0)
                return true;
        if(mpfr_cmp_si(a.right_mpfr(),b)<0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator==(const Gmpfi &a,int b){
        if(a.is_point())
                return(mpfr_cmp_si(a.left_mpfr(),b)?false:true);
        if(mpfr_cmp_si(a.left_mpfr(),b)>0 || mpfr_cmp_si(a.right_mpfr(),b)<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator<(const Gmpfi &a,long double b){
        if(mpfr_cmp_ld(a.right_mpfr(),b)<0)
                return true;
        if(mpfr_cmp_ld(a.left_mpfr(),b)>0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator>(const Gmpfi &a,long double b){
        if(mpfr_cmp_ld(a.left_mpfr(),b)>0)
                return true;
        if(mpfr_cmp_ld(a.right_mpfr(),b)<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator==(const Gmpfi &a,long double b){
        if(a.is_point())
                return(mpfr_cmp_ld(a.left_mpfr(),b)?false:true);
        if(mpfr_cmp_ld(a.left_mpfr(),b)>0 || mpfr_cmp_ld(a.right_mpfr(),b)<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator<(const Gmpfi &a,double b){
        if(mpfr_cmp_d(a.right_mpfr(),b)<0)
                return true;
        if(mpfr_cmp_d(a.left_mpfr(),b)>0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator>(const Gmpfi &a,double b){
        if(mpfr_cmp_d(a.left_mpfr(),b)>0)
                return true;
        if(mpfr_cmp_d(a.right_mpfr(),b)<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator==(const Gmpfi &a,double b){
        if(a.is_point())
                return(mpfr_cmp_d(a.left_mpfr(),b)?false:true);
        if(mpfr_cmp_d(a.left_mpfr(),b)>0 || mpfr_cmp_d(a.right_mpfr(),b)<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator<(const Gmpfi &a,const Gmpz &b){
        if(mpfr_cmp_z(a.right_mpfr(),b.mpz())<0)
                return true;
        if(mpfr_cmp_z(a.left_mpfr(),b.mpz())>0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator>(const Gmpfi &a,const Gmpz &b){
        if(mpfr_cmp_z(a.left_mpfr(),b.mpz())>0)
                return true;
        if(mpfr_cmp_z(a.right_mpfr(),b.mpz())<0 || a.is_point())
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<bool> operator==(const Gmpfi &a,const Gmpz &b){
        if(a.is_point())
                return(mpfr_cmp_z(a.left_mpfr(),b.mpz())?false:true);
        if(mpfr_cmp_z(a.left_mpfr(),b.mpz())>0 ||
                        mpfr_cmp_z(a.right_mpfr(),b.mpz())<0)
                return false;
        return Uncertain<bool>::indeterminate();
}

inline
Uncertain<Gmpfi> min BOOST_PREVENT_MACRO_SUBSTITUTION
                (const Gmpfi &x,const Gmpfi &y){
        return (x<=y)?x:y;
}

inline
Uncertain<Gmpfi> max BOOST_PREVENT_MACRO_SUBSTITUTION
                (const Gmpfi &x,const Gmpfi &y){
        return (x>=y)?x:y;
}

} // namespace CGAL

#endif  // CGAL_GMPFI_TYPE_H

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab

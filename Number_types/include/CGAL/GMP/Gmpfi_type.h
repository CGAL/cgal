// Copyright (c) 2007-2010 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_GMPFI_TYPE_H
#define CGAL_GMPFI_TYPE_H

#include <CGAL/config.h>
#include <CGAL/gmp.h>
#include <mpfr.h>
#include <CGAL/GMP/Gmpfr_type.h>
#include <CGAL/GMP/Gmpq_type.h>
#include <mpfi.h>
#include <boost/operators.hpp>
#include <CGAL/Uncertain.h>
#include <CGAL/tss.h>
#include <CGAL/IO/io.h>

#include <limits>
#include <algorithm>

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

// the default precision of Gmpfi is the size of a double's mantissa
#ifdef IEEE_DBL_MANT_DIG
#  define CGAL_GMPFI_DEFAULT_PRECISION IEEE_DBL_MANT_DIG
#else
#  define CGAL_GMPFI_DEFAULT_PRECISION 53
#endif

// the default precision is a variable local to each thread in multithreaded
// environments, or a global variable otherwise



class Gmpfi:
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
        private:

        // The endpoints of the interval are represented by two objects of
        // type Gmpfr. To apply MPFI functions to this interval, the
        // pointers to the data in _left and _right are copied to the
        // _interval structure using the function mpfi(). After the
        // operation, the function gather_bounds should be called to put
        // back the result of the operation in _left and _right.
        Gmpfr _left,_right;
        mutable __mpfi_struct _interval;

  static mp_prec_t&  default_precision()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(mp_prec_t, Gmpfi_default_precision, CGAL_GMPFI_DEFAULT_PRECISION);
    return Gmpfi_default_precision;
  }
        

        bool is_unique(){
#ifdef CGAL_GMPFR_NO_REFCOUNT
                return true;
#else
                return(_left.is_unique()&&_right.is_unique());
#endif
        }

        // swaps the contents of this object and another one
        void swap(Gmpfi &fi){
                std::swap(*this,fi);
        }

        // after calling a library function that modifies the data in the
        // structure _interval, this function has to be called in order to
        // copy the data in _interval to _left and _right
        void gather_bounds(){
                mpfr_custom_init_set(
                        _left.fr(),
                        mpfr_custom_get_kind(&_interval.left),
                        mpfr_custom_get_exp(&_interval.left),
                        mpfr_get_prec(&_interval.left),
                        mpfr_custom_get_mantissa(&_interval.left));
                mpfr_custom_init_set(
                        _right.fr(),
                        mpfr_custom_get_kind(&_interval.right),
                        mpfr_custom_get_exp(&_interval.right),
                        mpfr_get_prec(&_interval.right),
                        mpfr_custom_get_mantissa(&_interval.right));
        }

        public:

        typedef Gmpfr::Precision_type   Precision_type;

        // access

        mpfi_srcptr mpfi()const{
                _interval.left=*_left.fr();
                _interval.right=*_right.fr();
                CGAL_assertion(mpfr_equal_p(_left.fr(),&_interval.left)!=0 &&
                               mpfr_equal_p(_right.fr(),&_interval.right)!=0);
                return &_interval;
        }

        mpfi_ptr mpfi(){
                _interval.left=*_left.fr();
                _interval.right=*_right.fr();
                CGAL_assertion(mpfr_equal_p(_left.fr(),&_interval.left)!=0 &&
                               mpfr_equal_p(_right.fr(),&_interval.right)!=0);
                return &_interval;
        }

        mpfr_srcptr left_mpfr()const{
                return _left.fr();
        }

        mpfr_srcptr right_mpfr()const{
                return _right.fr();
        }

        Gmpfr inf()const{
                return _left;
        }

        Gmpfr sup()const{
                return _right;
        }

        // construction

        Gmpfi(){}
        ~Gmpfi(){}

#define CGAL_GMPFI_CONSTRUCTOR_FROM_SCALAR(_type) \
        Gmpfi(const _type &t, \
              Gmpfi::Precision_type p=Gmpfi::get_default_precision()){ \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                _left=Gmpfr(t,std::round_toward_neg_infinity,p); \
                _right=Gmpfr(t,std::round_toward_infinity,p); \
                CGAL_assertion(_left<=t&&_right>=t); \
        }

CGAL_GMPFI_CONSTRUCTOR_FROM_SCALAR(long);
CGAL_GMPFI_CONSTRUCTOR_FROM_SCALAR(unsigned long);
CGAL_GMPFI_CONSTRUCTOR_FROM_SCALAR(int);
CGAL_GMPFI_CONSTRUCTOR_FROM_SCALAR(double);
CGAL_GMPFI_CONSTRUCTOR_FROM_SCALAR(long double);
CGAL_GMPFI_CONSTRUCTOR_FROM_SCALAR(Gmpz);

#undef CGAL_GMPFI_CONSTRUCTOR_FROM_SCALAR

        Gmpfi(const Gmpq &q,
              Gmpfi::Precision_type p=Gmpfi::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                _left=Gmpfr(0,p);
                _right=Gmpfr(0,p);
                mpfr_set_q(_left.fr(),q.mpq(),GMP_RNDD);
                mpfr_set_q(_right.fr(),q.mpq(),GMP_RNDU);
                CGAL_assertion(_left<=q&&_right>=q);
        }

        Gmpfi(mpfi_srcptr i){
                _left=Gmpfr(&(i->left));
                _right=Gmpfr(&(i->right));
        }

        Gmpfi(mpfi_srcptr i,Gmpfi::Precision_type p){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                _left=Gmpfr(&(i->left),std::round_toward_neg_infinity,p);
                _right=Gmpfr(&(i->right),std::round_toward_infinity,p);
                CGAL_assertion(mpfr_cmp(_left.fr(),&(i->left))<=0 &&
                               mpfr_cmp(_right.fr(),&(i->right))>=0);
        }

        Gmpfi(const Gmpfr &f,
              Gmpfi::Precision_type p=Gmpfi::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                _left=Gmpfr(f,std::round_toward_neg_infinity,p);
                _right=Gmpfr(f,std::round_toward_infinity,p);
                CGAL_assertion(_left<=f&&_right>=f);
        }

        Gmpfi(const Gmpfr &l,
              const Gmpfr &r,
              Gmpfi::Precision_type p=Gmpfi::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                _left=Gmpfr(l,std::round_toward_neg_infinity,p);
                _right=Gmpfr(r,std::round_toward_infinity,p);
                CGAL_assertion(_left<=l||(_left.is_nan()&&l.is_nan()));
                CGAL_assertion(_right>=l||(_right.is_nan()&&r.is_nan()));
        }

        Gmpfi(const std::pair<Gmpfr,Gmpfr> &bounds,
              Gmpfi::Precision_type p=Gmpfi::get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                _left=Gmpfr(bounds.first,std::round_toward_neg_infinity,p);
                _right=Gmpfr(bounds.second,std::round_toward_infinity,p);
                CGAL_assertion(_left<=bounds.first||
                               (_left.is_nan()&&bounds.first.is_nan()));
                CGAL_assertion(_right>=bounds.second||
                               (_right.is_nan()&&bounds.second.is_nan()));
        }

        template<class L,class R>
        Gmpfi(const std::pair<L,R> &bounds,
              Gmpfi::Precision_type p=get_default_precision()){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                _left=Gmpfr(bounds.first,std::round_toward_neg_infinity,p);
                _right=Gmpfr(bounds.second,std::round_toward_infinity,p);
                CGAL_assertion(_left<=bounds.first&&_right>=bounds.second);
        }

        // copy assignment operator
        Gmpfi& operator=(const Gmpfi &a){
                _left=a.inf();
                _right=a.sup();
                CGAL_assertion(_left==a.inf()||
                               (_left.is_nan()&&a.inf().is_nan()));
                CGAL_assertion(_right==a.sup()||
                               (_right.is_nan()&&a.sup().is_nan()));
                return *this;
        }

        // copy constructor without precision
        Gmpfi(const Gmpfi &a){
                _left=a.inf();
                _right=a.sup();
                CGAL_assertion(_left==a.inf()||
                               (_left.is_nan()&&a.inf().is_nan()));
                CGAL_assertion(_right==a.sup()||
                               (_right.is_nan()&&a.sup().is_nan()));
        }

        // copy constructor with precision
        Gmpfi(const Gmpfi &a,Gmpfi::Precision_type p){
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
                _left=Gmpfr(a.inf(),std::round_toward_neg_infinity,p);
                _right=Gmpfr(a.sup(),std::round_toward_infinity,p);
                CGAL_assertion(_left<=a.inf()||
                               (_left.is_nan()&&a.inf().is_nan()));
                CGAL_assertion(_right>=a.sup()||
                               (_right.is_nan()&&a.sup().is_nan()));
        }

        // default precision


        static Gmpfi::Precision_type get_default_precision();
        static Gmpfi::Precision_type set_default_precision(
                                                Gmpfi::Precision_type prec);

        // precision of a single Gmpfi object

        Gmpfi::Precision_type get_precision()const;
        Gmpfi round(Gmpfi::Precision_type)const;

        // arithmetics

        Gmpfi operator+()const;
        Gmpfi operator-()const;

#define CGAL_GMPFI_DECLARE_OPERATORS(_type) \
        Gmpfi& operator+=(_type); \
        Gmpfi& operator-=(_type); \
        Gmpfi& operator*=(_type); \
        Gmpfi& operator/=(_type);

        CGAL_GMPFI_DECLARE_OPERATORS(const Gmpfi&)
        CGAL_GMPFI_DECLARE_OPERATORS(const Gmpfr&)
        CGAL_GMPFI_DECLARE_OPERATORS(long)
        CGAL_GMPFI_DECLARE_OPERATORS(unsigned long)
        CGAL_GMPFI_DECLARE_OPERATORS(int)
        CGAL_GMPFI_DECLARE_OPERATORS(const Gmpz&)
        CGAL_GMPFI_DECLARE_OPERATORS(const Gmpq&)

#undef CGAL_GMPFI_DECLARE_OPERATORS

#define CGAL_GMPFI_DECLARE_STATIC_FUNCTION(_f,_t1,_t2) \
        static Gmpfi _f (_t1,_t2,Gmpfi::Precision_type=0);

#define CGAL_GMPFI_DECLARE_STATIC_FUNCTIONS(_t1,_t2) \
        CGAL_GMPFI_DECLARE_STATIC_FUNCTION(add,_t1,_t2) \
        CGAL_GMPFI_DECLARE_STATIC_FUNCTION(sub,_t1,_t2) \
        CGAL_GMPFI_DECLARE_STATIC_FUNCTION(mul,_t1,_t2) \
        CGAL_GMPFI_DECLARE_STATIC_FUNCTION(div,_t1,_t2)

#define CGAL_GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(_t) \
        CGAL_GMPFI_DECLARE_STATIC_FUNCTIONS(const Gmpfi&,_t) \
        CGAL_GMPFI_DECLARE_STATIC_FUNCTIONS(_t,const Gmpfi&)

        CGAL_GMPFI_DECLARE_STATIC_FUNCTIONS(const Gmpfi&,const Gmpfi&)
        CGAL_GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(const Gmpfr&)
        CGAL_GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(long)
        CGAL_GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(unsigned long)
        CGAL_GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(int)
        CGAL_GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(const Gmpz&)
        CGAL_GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS(const Gmpq&)

#undef CGAL_GMPFI_DECLARE_STATIC_FUNCTION
#undef CGAL_GMPFI_DECLARE_STATIC_FUNCTIONS
#undef CGAL_GMPFI_DECLARE_TWO_WAY_STATIC_FUNCTIONS

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


inline
Gmpfi::Precision_type Gmpfi::get_default_precision(){

  return default_precision();
}

inline
Gmpfi::Precision_type Gmpfi::set_default_precision(Gmpfi::Precision_type prec){
        Gmpfi::Precision_type old_prec= default_precision();
        CGAL_assertion(prec>=MPFR_PREC_MIN&&prec<=MPFR_PREC_MAX);
        default_precision() = prec;

        return old_prec;
}

// precision of a single Gmpfi object

inline
Gmpfi::Precision_type Gmpfi::get_precision()const{
        return (_left.get_precision()>_right.get_precision()?
                (Gmpfi::Precision_type)_left.get_precision():
                (Gmpfi::Precision_type)_right.get_precision());
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
        Gmpfi result (0, this->get_precision());
        mpfi_neg (result.mpfi(), this->mpfi());
        *(result._left.fr()) = result._interval.left;
        *(result._right.fr()) = result._interval.right;
        return result;
}

// CGAL_GMPFI_BALANCE_ENDPOINTS checks if both bounds of the interval have
// the same precision. If not, it rounds the one with the smallest
// precision.
#define CGAL_GMPFI_BALANCE_ENDPOINTS \
        if(_left.get_precision()<_right.get_precision()){ \
                _left=Gmpfr(_left,_right.get_precision()); \
        }else{ \
                if(_right.get_precision()<_left.get_precision()){ \
                        _right=Gmpfr(_right,_left.get_precision()); \
                } \
        }; \
        CGAL_assertion_msg(_left.get_precision()==_right.get_precision(), \
                           "error balancing bounds precision");

// CGAL_GMPFI_OBJECT_BINARY_OPERATOR defines an overloaded binary operator
// of the Gmpfi class, where the operated object belongs to another class,
// which represents a point (as opposition to an interval). The operation
// will be performed using the biggest precision of the bounds of this
// Gmpfi object. That means that if bounds have different precision, one
// of them (the one with the smallest precision) will be rounded. This is
// not a problem when the object is not unique, since a new Gmpfi object
// will be created with the bounds having the correct precision.
#define CGAL_GMPFI_OBJECT_BINARY_OPERATOR(_op,_class,_member,_fun) \
        inline \
        Gmpfi& Gmpfi::_op(const _class &b){ \
                if(is_unique()){ \
                        CGAL_GMPFI_BALANCE_ENDPOINTS \
                        _fun(mpfi(),mpfi(),b._member); \
                        gather_bounds(); \
                }else{ \
                        Gmpfi result (0, this->get_precision()); \
                        _fun (result.mpfi(), this->mpfi(), b._member); \
                        *(result._left.fr()) = result._interval.left; \
                        *(result._right.fr()) = result._interval.right; \
                        this->swap (result); \
                } \
                return(*this); \
        }

// CGAL_GMPFI_GMPFI_BINARY_OPERATOR defines an overloaded binary operator
// of the Gmpfi class, where both operands are Gmpfi objects.  The
// operation will be performed using the biggest precision of the bounds
// of both intervals. The bounds of target object will be rounded
// accordingly before the operation.
#define CGAL_GMPFI_GMPFI_BINARY_OPERATOR(_op,_fun) \
        inline \
        Gmpfi& Gmpfi::_op(const Gmpfi &fi){ \
                if(is_unique()){ \
                        if(get_precision()<fi.get_precision()){ \
                                Gmpfi result (0, fi.get_precision()); \
                                _fun(result.mpfi(), this->mpfi(), fi.mpfi()); \
                                *(result._left.fr()) = result._interval.left; \
                                *(result._right.fr())= result._interval.right;\
                                this->swap (result); \
                        }else{ \
                                CGAL_GMPFI_BALANCE_ENDPOINTS \
                                _fun(mpfi(),mpfi(),fi.mpfi()); \
                                gather_bounds(); \
                        } \
                }else{ \
                        Gmpfi result(0, \
                                     this->get_precision()<fi.get_precision()?\
                                     fi.get_precision(): \
                                     this->get_precision()); \
                        _fun (result.mpfi(), this->mpfi(), fi.mpfi()); \
                        *(result._left.fr()) = result._interval.left; \
                        *(result._right.fr()) = result._interval.right; \
                        this->swap (result); \
                } \
                return(*this); \
  }

// CGAL_GMPFI_TYPE_BINARY_OPERATOR defines an overloaded binary operator of
// the Gmpfi class, where the operated number belongs to a c++ type.
// Precision of the operation is defined in the same manner that in
// CGAL_GMPFI_OBJECT_BINARY_OPERATOR.
#define CGAL_GMPFI_TYPE_BINARY_OPERATOR(_op,_type,_fun) \
        inline \
        Gmpfi& Gmpfi::_op(_type x){ \
                if(is_unique()){ \
                        CGAL_GMPFI_BALANCE_ENDPOINTS \
                        _fun(mpfi(),mpfi(),x); \
                        gather_bounds(); \
                }else{ \
                        Gmpfi result (0, this->get_precision()); \
                        _fun (result.mpfi(), this->mpfi(), x); \
                        *(result._left.fr()) = result._interval.left; \
                        *(result._right.fr()) = result._interval.right; \
                        this->swap (result); \
                } \
                return *this; \
        }

CGAL_GMPFI_GMPFI_BINARY_OPERATOR(operator+=,mpfi_add)
CGAL_GMPFI_GMPFI_BINARY_OPERATOR(operator-=,mpfi_sub)
CGAL_GMPFI_GMPFI_BINARY_OPERATOR(operator*=,mpfi_mul)
CGAL_GMPFI_GMPFI_BINARY_OPERATOR(operator/=,mpfi_div)

CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator+=,Gmpfr,fr(),mpfi_add_fr)
CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator-=,Gmpfr,fr(),mpfi_sub_fr)
CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator*=,Gmpfr,fr(),mpfi_mul_fr)
CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator/=,Gmpfr,fr(),mpfi_div_fr)

CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator+=,long,mpfi_add_si)
CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator-=,long,mpfi_sub_si)
CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator*=,long,mpfi_mul_si)
CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator/=,long,mpfi_div_si)

CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator+=,unsigned long,mpfi_add_ui)
CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator-=,unsigned long,mpfi_sub_ui)
CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator*=,unsigned long,mpfi_mul_ui)
CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator/=,unsigned long,mpfi_div_ui)

CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator+=,int,mpfi_add_si)
CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator-=,int,mpfi_sub_si)
CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator*=,int,mpfi_mul_si)
CGAL_GMPFI_TYPE_BINARY_OPERATOR(operator/=,int,mpfi_div_si)

CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator+=,Gmpz,mpz(),mpfi_add_z)
CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator-=,Gmpz,mpz(),mpfi_sub_z)
CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator*=,Gmpz,mpz(),mpfi_mul_z)
CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator/=,Gmpz,mpz(),mpfi_div_z)

CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator+=,Gmpq,mpq(),mpfi_add_q)
CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator-=,Gmpq,mpq(),mpfi_sub_q)
CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator*=,Gmpq,mpq(),mpfi_mul_q)
CGAL_GMPFI_OBJECT_BINARY_OPERATOR(operator/=,Gmpq,mpq(),mpfi_div_q)

#undef CGAL_GMPFI_BALANCE_ENDPOINTS
#undef CGAL_GMPFI_GMPFI_BINARY_OPERATOR
#undef CGAL_GMPFI_OBJECT_BINARY_OPERATOR
#undef CGAL_GMPFI_TYPE_BINARY_OPERATOR

// the static arithmetic functions are defined in a separate file
#include <CGAL/GMP/Gmpfi_type_static.h>

#define CGAL_GMPFI_ARITHMETIC_FUNCTION(_name,_fun) \
        inline \
        Gmpfi Gmpfi::_name (Gmpfi::Precision_type p)const{ \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result (0, p); \
                _fun (result.mpfi(), this->mpfi()); \
                *(result._left.fr()) = result._interval.left; \
                *(result._right.fr()) = result._interval.right; \
                return result; \
        }

CGAL_GMPFI_ARITHMETIC_FUNCTION(abs,mpfi_abs)
CGAL_GMPFI_ARITHMETIC_FUNCTION(sqrt,mpfi_sqrt)

inline
Gmpfi Gmpfi::cbrt(Gmpfi::Precision_type p)const{
        // MPFI does not provide a cubic root function
        CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
        Gmpfi result (0, p);
        mpfr_cbrt(&(result.mpfi())->left, left_mpfr(), GMP_RNDD);
        mpfr_cbrt(&(result.mpfi())->right,right_mpfr(),GMP_RNDU);
        *(result._left.fr()) = result._interval.left;
        *(result._right.fr()) = result._interval.right;
        return result;
}

inline
Gmpfi Gmpfi::kthroot(int k,Gmpfi::Precision_type p)const{
        // MPFI does not provide k-th root functions
        CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX);
        Gmpfi result (0, p);
        mpfr_root(&(result.mpfi())->left, left_mpfr(), k,GMP_RNDD);
        mpfr_root(&(result.mpfi())->right,right_mpfr(),k,GMP_RNDU);
        *(result._left.fr()) = result._interval.left;
        *(result._right.fr()) = result._interval.right;
        return result;
}

CGAL_GMPFI_ARITHMETIC_FUNCTION(square,mpfi_sqr)

#undef CGAL_GMPFI_ARITHMETIC_FUNCTION

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
        int leftsign=mpfr_sgn(left_mpfr());
        if(leftsign>0)
                return POSITIVE;
        if(leftsign==0){
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
        if(mpfr_zero_p(left_mpfr())!=0 && mpfr_zero_p(right_mpfr())!=0)
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
        // if a bound is finite and its double is infinity, we overflow
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
        long *e=NULL;
        mpfr_init2(middle,53);
        mpfi_get_fr(middle,mpfi());
        double d=mpfr_get_d_2exp(e,middle,mpfr_get_default_rounding_mode());
        mpfr_clear(middle);
        return std::make_pair(d,*e);
}

inline
std::pair<std::pair<double,double>,long> Gmpfi::to_interval_exp()const{
        long *e1=NULL,*e2=NULL;
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
        internal::eat_white_space(is);
        c=is.get();
        if(c!='['){
                invalid_number:
                is.setstate(std::ios_base::failbit);
                is.flags(old_flags);
                return is;
        }
        internal::eat_white_space(is);
        is>>left;
        c=is.get();
        if(c!=',')
                goto invalid_number;
        is>>right;
        internal::eat_white_space(is);
        c=is.get();
        if(c!=']')
                goto invalid_number;
        // Why is this done the following way? Because left and right can
        // have different precision. Doing this with a constructor would
        // force to create a Gmpfi where both endpoints have the same
        // precision, what can give a wrong reconstruction of a previously
        // outputted number. (This function will give a good reconstruction
        // iff Gmpfr gives a good reconstruction.)
        Gmpfi temp(0,(Gmpfi::Precision_type)MPFR_PREC_MIN);
        mpfr_swap(left.fr(), temp.inf().fr());
        mpfr_swap(right.fr(),temp.sup().fr());
        f=temp;
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

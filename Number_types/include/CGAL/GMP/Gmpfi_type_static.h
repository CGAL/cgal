// Copyright (c) 2007-2010 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

// This file contains the arithmetic functions not members of the Gmpfi
// class.

// CGAL_GMPFI_PREC returns the precision used to operate a Gmpfi object _a
// and an object of another type. Currently, the returned value is the
// precision of _a.
#define CGAL_GMPFI_PREC(_a) ( _a.get_precision() )

// CGAL_GMPFI_PREC_2 returns the precision used to operate between two
// Gmpfi objects _a and _b. Currently, the returned value is the maximum
// between the precisions of _a and _b.
#define CGAL_GMPFI_PREC_2(_a,_b) \
        ( _a.get_precision() > _b.get_precision() ? \
          _a.get_precision() : \
          _b.get_precision() )

// CGAL_GMPFI_OP_GMPFI defines an arithmetic operation between two
// Gmpfi objects.
#define CGAL_GMPFI_OP_GMPFI(_name,_fun) \
        inline \
        Gmpfi Gmpfi::_name (const Gmpfi &a, \
                            const Gmpfi &b, \
                            Gmpfi::Precision_type p){ \
                if(!p) \
                        p=CGAL_GMPFI_PREC_2(a,b); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result (0,p); \
                _fun(result.mpfi(),a.mpfi(),b.mpfi()); \
                *(result._left.fr())=result._interval.left; \
                *(result._right.fr())=result._interval.right; \
                return result; \
}

// CGAL_GMPFI_COMMUTATIVE_OP defines a commutative arithmetic operation
// between a Gmpfi object and an object of another type.
#define CGAL_GMPFI_COMMUTATIVE_OP(_name,_type,_member,_fun) \
        inline \
        Gmpfi Gmpfi::_name (const Gmpfi &a,_type b,Gmpfi::Precision_type p){ \
                if(!p) \
                        p=CGAL_GMPFI_PREC(a); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result(0,p); \
                _fun(result.mpfi(),a.mpfi(),_member); \
                *(result._left.fr())=result._interval.left; \
                *(result._right.fr())=result._interval.right; \
                return result; \
        } \
        inline \
        Gmpfi Gmpfi::_name (_type b,const Gmpfi &a,Gmpfi::Precision_type p){ \
                if(!p) \
                        p=CGAL_GMPFI_PREC(a); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result(0,p); \
                _fun(result.mpfi(),a.mpfi(),_member); \
                *(result._left.fr())=result._interval.left; \
                *(result._right.fr())=result._interval.right; \
                return result; \
        }

// CGAL_GMPFI_NONCOMMUTATIVE_OP defines a non-commutative arithmetic
// operation between a Gmpfi object and an object of another type.
#define CGAL_GMPFI_NONCOMMUTATIVE_OP(_name,_type,_member,_fun1,_fun2) \
        inline \
        Gmpfi Gmpfi::_name (const Gmpfi &a,_type b,Gmpfi::Precision_type p){ \
                if(!p) \
                        p=CGAL_GMPFI_PREC(a); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result(0,p); \
                _fun1(result.mpfi(),a.mpfi(),_member); \
                *(result._left.fr())=result._interval.left; \
                *(result._right.fr())=result._interval.right; \
                return result; \
        } \
        inline \
        Gmpfi Gmpfi::_name (_type b,const Gmpfi &a,Gmpfi::Precision_type p){ \
                if(!p) \
                        p=CGAL_GMPFI_PREC(a); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result(0,p); \
                _fun2(result.mpfi(),_member,a.mpfi()); \
                *(result._left.fr())=result._interval.left; \
                *(result._right.fr())=result._interval.right; \
                return result; \
        }

CGAL_GMPFI_OP_GMPFI(add,mpfi_add)
CGAL_GMPFI_OP_GMPFI(sub,mpfi_sub)
CGAL_GMPFI_OP_GMPFI(mul,mpfi_mul)
CGAL_GMPFI_OP_GMPFI(div,mpfi_div)

CGAL_GMPFI_COMMUTATIVE_OP(add,const Gmpfr&,b.fr(),mpfi_add_fr)
CGAL_GMPFI_NONCOMMUTATIVE_OP(sub,const Gmpfr&,b.fr(),mpfi_sub_fr,mpfi_fr_sub)
CGAL_GMPFI_COMMUTATIVE_OP(mul,const Gmpfr&,b.fr(),mpfi_mul_fr)
CGAL_GMPFI_NONCOMMUTATIVE_OP(div,const Gmpfr&,b.fr(),mpfi_div_fr,mpfi_fr_div)

CGAL_GMPFI_COMMUTATIVE_OP(add,long,b,mpfi_add_si)
CGAL_GMPFI_NONCOMMUTATIVE_OP(sub,long,b,mpfi_sub_si,mpfi_si_sub)
CGAL_GMPFI_COMMUTATIVE_OP(mul,long,b,mpfi_mul_si)
CGAL_GMPFI_NONCOMMUTATIVE_OP(div,long,b,mpfi_div_si,mpfi_si_div)

CGAL_GMPFI_COMMUTATIVE_OP(add,unsigned long,b,mpfi_add_ui)
CGAL_GMPFI_NONCOMMUTATIVE_OP(sub,unsigned long,b,mpfi_sub_ui,mpfi_ui_sub)
CGAL_GMPFI_COMMUTATIVE_OP(mul,unsigned long,b,mpfi_mul_ui)
CGAL_GMPFI_NONCOMMUTATIVE_OP(div,unsigned long,b,mpfi_div_ui,mpfi_ui_div)

CGAL_GMPFI_COMMUTATIVE_OP(add,int,b,mpfi_add_si)
CGAL_GMPFI_NONCOMMUTATIVE_OP(sub,int,b,mpfi_sub_si,mpfi_si_sub)
CGAL_GMPFI_COMMUTATIVE_OP(mul,int,b,mpfi_mul_si)
CGAL_GMPFI_NONCOMMUTATIVE_OP(div,int,b,mpfi_div_si,mpfi_si_div)

CGAL_GMPFI_COMMUTATIVE_OP(add,const Gmpz&,b.mpz(),mpfi_add_z)
CGAL_GMPFI_NONCOMMUTATIVE_OP(sub,const Gmpz&,b.mpz(),mpfi_sub_z,mpfi_z_sub)
CGAL_GMPFI_COMMUTATIVE_OP(mul,const Gmpz&,b.mpz(),mpfi_mul_z)
CGAL_GMPFI_NONCOMMUTATIVE_OP(div,const Gmpz&,b.mpz(),mpfi_div_z,mpfi_z_div)

CGAL_GMPFI_COMMUTATIVE_OP(add,const Gmpq&,b.mpq(),mpfi_add_q)
CGAL_GMPFI_NONCOMMUTATIVE_OP(sub,const Gmpq&,b.mpq(),mpfi_sub_q,mpfi_q_sub)
CGAL_GMPFI_COMMUTATIVE_OP(mul,const Gmpq&,b.mpq(),mpfi_mul_q)
CGAL_GMPFI_NONCOMMUTATIVE_OP(div,const Gmpq&,b.mpq(),mpfi_div_q,mpfi_q_div)

#undef CGAL_GMPFI_PREC
#undef CGAL_GMPFI_PREC_2
#undef CGAL_GMPFI_OP_GMPFI
#undef CGAL_GMPFI_COMMUTATIVE_OP
#undef CGAL_GMPFI_NONCOMMUTATIVE_OP

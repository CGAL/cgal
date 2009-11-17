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

// This file contains the arithmetic functions not members of the Gmpfi
// class.

// _GMPFI_PREC returns the precision used to operate a Gmpfi object _a
// and an object of another type. Currently, the returned value is the
// precision of _a.
#define _GMPFI_PREC(_a) ( _a.get_precision() )

// _GMPFI_PREC_2 returns the precision used to operate between two
// Gmpfi objects _a and _b. Currently, the returned value is the maximum
// between the precisions of _a and _b.
#define _GMPFI_PREC_2(_a,_b) \
        ( _a.get_precision() > _b.get_precision() ? \
          _a.get_precision() : \
          _b.get_precision() )

// _GMPFI_OPERATION_GMPFI defines an arithmetic operation between two
// Gmpfi objects.
#define _GMPFI_OPERATION_GMPFI(_name,_fun) \
        inline \
        Gmpfi Gmpfi::_name (const Gmpfi &a, \
                            const Gmpfi &b, \
                            Gmpfi::Precision_type p){ \
                if(!p) \
                        p=_GMPFI_PREC_2(a,b); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result(0,p); \
                _fun(result.mpfi(),a.mpfi(),b.mpfi()); \
                return result; \
}

// _GMPFI_COMMUTATIVE_OPERATION defines a commutative arithmetic operation
// between a Gmpfi object and an object of another type.
#define _GMPFI_COMMUTATIVE_OPERATION(_name,_type,_member,_fun) \
        inline \
        Gmpfi Gmpfi::_name (const Gmpfi &a,_type b,Gmpfi::Precision_type p){ \
                if(!p) \
                        p=_GMPFI_PREC(a); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result(0,p); \
                _fun(result.mpfi(),a.mpfi(),_member); \
                return result; \
        } \
        inline \
        Gmpfi Gmpfi::_name (_type b,const Gmpfi &a,Gmpfi::Precision_type p){ \
                if(!p) \
                        p=_GMPFI_PREC(a); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result(0,p); \
                _fun(result.mpfi(),a.mpfi(),_member); \
                return result; \
        }

// _GMPFI_NONCOMMUTATIVE_OPERATION defines a non-commutative arithmetic
// operation between a Gmpfi object and an object of another type.
#define _GMPFI_NONCOMMUTATIVE_OPERATION(_name,_type,_member,_fun1,_fun2) \
        inline \
        Gmpfi Gmpfi::_name (const Gmpfi &a,_type b,Gmpfi::Precision_type p){ \
                if(!p) \
                        p=_GMPFI_PREC(a); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result(0,p); \
                _fun1(result.mpfi(),a.mpfi(),_member); \
                return result; \
        } \
        inline \
        Gmpfi Gmpfi::_name (_type b,const Gmpfi &a,Gmpfi::Precision_type p){ \
                if(!p) \
                        p=_GMPFI_PREC(a); \
                CGAL_assertion(p>=MPFR_PREC_MIN&&p<=MPFR_PREC_MAX); \
                Gmpfi result(0,p); \
                _fun2(result.mpfi(),_member,a.mpfi()); \
                return result; \
        }

_GMPFI_OPERATION_GMPFI(add,mpfi_add)
_GMPFI_OPERATION_GMPFI(sub,mpfi_sub)
_GMPFI_OPERATION_GMPFI(mul,mpfi_mul)
_GMPFI_OPERATION_GMPFI(div,mpfi_div)

_GMPFI_COMMUTATIVE_OPERATION(add,const Gmpfr&,b.fr(),mpfi_add_fr)
_GMPFI_NONCOMMUTATIVE_OPERATION(sub,const Gmpfr&,b.fr(),mpfi_sub_fr,mpfi_fr_sub)
_GMPFI_COMMUTATIVE_OPERATION(mul,const Gmpfr&,b.fr(),mpfi_mul_fr)
_GMPFI_NONCOMMUTATIVE_OPERATION(div,const Gmpfr&,b.fr(),mpfi_div_fr,mpfi_fr_div)

_GMPFI_COMMUTATIVE_OPERATION(add,long,b,mpfi_add_si)
_GMPFI_NONCOMMUTATIVE_OPERATION(sub,long,b,mpfi_sub_si,mpfi_si_sub)
_GMPFI_COMMUTATIVE_OPERATION(mul,long,b,mpfi_mul_si)
_GMPFI_NONCOMMUTATIVE_OPERATION(div,long,b,mpfi_div_si,mpfi_si_div)

_GMPFI_COMMUTATIVE_OPERATION(add,unsigned long,b,mpfi_add_ui)
_GMPFI_NONCOMMUTATIVE_OPERATION(sub,unsigned long,b,mpfi_sub_ui,mpfi_ui_sub)
_GMPFI_COMMUTATIVE_OPERATION(mul,unsigned long,b,mpfi_mul_ui)
_GMPFI_NONCOMMUTATIVE_OPERATION(div,unsigned long,b,mpfi_div_ui,mpfi_ui_div)

_GMPFI_COMMUTATIVE_OPERATION(add,int,b,mpfi_add_si)
_GMPFI_NONCOMMUTATIVE_OPERATION(sub,int,b,mpfi_sub_si,mpfi_si_sub)
_GMPFI_COMMUTATIVE_OPERATION(mul,int,b,mpfi_mul_si)
_GMPFI_NONCOMMUTATIVE_OPERATION(div,int,b,mpfi_div_si,mpfi_si_div)

_GMPFI_COMMUTATIVE_OPERATION(add,const Gmpz&,b.mpz(),mpfi_add_z)
_GMPFI_NONCOMMUTATIVE_OPERATION(sub,const Gmpz&,b.mpz(),mpfi_sub_z,mpfi_z_sub)
_GMPFI_COMMUTATIVE_OPERATION(mul,const Gmpz&,b.mpz(),mpfi_mul_z)
_GMPFI_NONCOMMUTATIVE_OPERATION(div,const Gmpz&,b.mpz(),mpfi_div_z,mpfi_z_div)

_GMPFI_COMMUTATIVE_OPERATION(add,const Gmpq&,b.mpq(),mpfi_add_q)
_GMPFI_NONCOMMUTATIVE_OPERATION(sub,const Gmpq&,b.mpq(),mpfi_sub_q,mpfi_q_sub)
_GMPFI_COMMUTATIVE_OPERATION(mul,const Gmpq&,b.mpq(),mpfi_mul_q)
_GMPFI_NONCOMMUTATIVE_OPERATION(div,const Gmpq&,b.mpq(),mpfi_div_q,mpfi_q_div)

#undef _GMPFI_PREC
#undef _GMPFI_PREC_2 
#undef _GMPFI_OPERATION_GMPFI
#undef _GMPFI_COMMUTATIVE_OPERATION
#undef _GMPFI_NONCOMMUTATIVE_OPERATION

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab

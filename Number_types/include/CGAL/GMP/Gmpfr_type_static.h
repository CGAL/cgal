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

// This file contains the arithmetic functions not members of the Gmpfr
// class.

// _GMPFR_PREC returns the precision used to operate a Gmpfr object _a
// and an object of another type. Currently, the returned value is the
// maximum between the precision of _a and the default precision.
#define _GMPFR_PREC(_a) \
        ( mpfr_get_prec(_a.fr()) > Gmpfr::get_default_precision() ? \
          mpfr_get_prec(_a.fr()): \
          Gmpfr::get_default_precision() )

// _GMPFR_PREC_2 returns the precision used to operate between two
// Gmpfr objects _a and _b. Currently, the returned value is the maximum
// between the precisions of _a and _b and the default precision. Above
// comment on optimality also holds for this function.
#define _GMPFR_PREC_2(_a,_b) \
        ( mpfr_get_prec(_a.fr()) >= mpfr_get_prec(_b.fr()) ? \
          ( mpfr_get_prec(_a.fr()) > Gmpfr::get_default_precision() ? \
            mpfr_get_prec(_a.fr()): \
            Gmpfr::get_default_precision() ) : \
          ( mpfr_get_prec(_b.fr()) > Gmpfr::get_default_precision() ? \
            mpfr_get_prec(_b.fr()): \
            Gmpfr::get_default_precision() ) )

// _GMPFR_OPERATION_GMPFR defines an arithmetic operation between two
// Gmpfr objects.
#define _GMPFR_OPERATION_GMPFR(_name,_fun) \
        inline \
        Gmpfr Gmpfr::_name (const Gmpfr &a,const Gmpfr &b,std::float_round_style r){ \
                Gmpfr result(0,_GMPFR_PREC_2(a,b)); \
                _fun(result.fr(),a.fr(),b.fr(),_gmp_rnd(r)); \
                return result; \
        } \
        inline \
        Gmpfr Gmpfr::_name (const Gmpfr &a, \
                            const Gmpfr &b, \
                            Gmpfr::Precision_type p, \
                            std::float_round_style r){ \
                Gmpfr result(0,p); \
                _fun(result.fr(),a.fr(),b.fr(),_gmp_rnd(r)); \
                return result; \
}

// _GMPFR_NONCOMMUTATIVE_OPERATION defines an arithmetic operation between
// a Gmpfr object and an object of another type. All operations are treated
// as non-commutative, since MPFR does not provide all the functions needed
// to operate with any type as first operand.
#define _GMPFR_NONCOMMUTATIVE_OPERATION(_name,_type,_member,_fun) \
        inline \
        Gmpfr Gmpfr::_name (const Gmpfr &a,_type b,std::float_round_style r){ \
                Gmpfr result(0,_GMPFR_PREC(a)); \
                _fun(result.fr(),a.fr(),_member,_gmp_rnd(r)); \
                return result; \
        } \
        inline \
        Gmpfr Gmpfr::_name (const Gmpfr &a, \
                            _type b, \
                            Gmpfr::Precision_type p, \
                            std::float_round_style r){ \
                Gmpfr result(0,p); \
                _fun(result.fr(),a.fr(),_member,_gmp_rnd(r)); \
                return result; \
        }

_GMPFR_OPERATION_GMPFR(add,mpfr_add)
_GMPFR_OPERATION_GMPFR(sub,mpfr_sub)
_GMPFR_OPERATION_GMPFR(mul,mpfr_mul)
_GMPFR_OPERATION_GMPFR(div,mpfr_div)

_GMPFR_NONCOMMUTATIVE_OPERATION(add,long,b,mpfr_add_si)
_GMPFR_NONCOMMUTATIVE_OPERATION(sub,long,b,mpfr_sub_si)
_GMPFR_NONCOMMUTATIVE_OPERATION(mul,long,b,mpfr_mul_si)
_GMPFR_NONCOMMUTATIVE_OPERATION(div,long,b,mpfr_div_si)

_GMPFR_NONCOMMUTATIVE_OPERATION(add,unsigned long,b,mpfr_add_ui)
_GMPFR_NONCOMMUTATIVE_OPERATION(sub,unsigned long,b,mpfr_sub_ui)
_GMPFR_NONCOMMUTATIVE_OPERATION(mul,unsigned long,b,mpfr_mul_ui)
_GMPFR_NONCOMMUTATIVE_OPERATION(div,unsigned long,b,mpfr_div_ui)

_GMPFR_NONCOMMUTATIVE_OPERATION(add,int,b,mpfr_add_si)
_GMPFR_NONCOMMUTATIVE_OPERATION(sub,int,b,mpfr_sub_si)
_GMPFR_NONCOMMUTATIVE_OPERATION(mul,int,b,mpfr_mul_si)
_GMPFR_NONCOMMUTATIVE_OPERATION(div,int,b,mpfr_div_si)

_GMPFR_NONCOMMUTATIVE_OPERATION(add,const Gmpz&,b.mpz(),mpfr_add_z)
_GMPFR_NONCOMMUTATIVE_OPERATION(sub,const Gmpz&,b.mpz(),mpfr_sub_z)
_GMPFR_NONCOMMUTATIVE_OPERATION(mul,const Gmpz&,b.mpz(),mpfr_mul_z)
_GMPFR_NONCOMMUTATIVE_OPERATION(div,const Gmpz&,b.mpz(),mpfr_div_z)

_GMPFR_NONCOMMUTATIVE_OPERATION(add,const Gmpq&,b.mpq(),mpfr_add_q)
_GMPFR_NONCOMMUTATIVE_OPERATION(sub,const Gmpq&,b.mpq(),mpfr_sub_q)
_GMPFR_NONCOMMUTATIVE_OPERATION(mul,const Gmpq&,b.mpq(),mpfr_mul_q)
_GMPFR_NONCOMMUTATIVE_OPERATION(div,const Gmpq&,b.mpq(),mpfr_div_q)

#undef _GMPFR_PREC
#undef _GMPFR_PREC_2 
#undef _GMPFR_OPERATION_GMPFR
#undef _GMPFR_NONCOMMUTATIVE_OPERATION

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab

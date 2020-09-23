// Copyright (c) 2007-2010 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

// This file contains the arithmetic functions not members of the Gmpfr
// class.

// CGAL_GMPFR_PREC returns the precision used to operate a Gmpfr object _a
// and an object of another type. Currently, the returned value is the
// maximum between the precision of _a and the default precision.
#define CGAL_GMPFR_PREC(_a) \
        ( mpfr_get_prec(_a.fr()) > Gmpfr::get_default_precision() ? \
          mpfr_get_prec(_a.fr()): \
          Gmpfr::get_default_precision() )

// CGAL_GMPFR_PREC_2 returns the precision used to operate between two
// Gmpfr objects _a and _b. Currently, the returned value is the maximum
// between the precisions of _a and _b and the default precision. Above
// comment on optimality also holds for this function.
#define CGAL_GMPFR_PREC_2(_a,_b) \
        ( mpfr_get_prec(_a.fr()) >= mpfr_get_prec(_b.fr()) ? \
          ( mpfr_get_prec(_a.fr()) > Gmpfr::get_default_precision() ? \
            mpfr_get_prec(_a.fr()): \
            Gmpfr::get_default_precision() ) : \
          ( mpfr_get_prec(_b.fr()) > Gmpfr::get_default_precision() ? \
            mpfr_get_prec(_b.fr()): \
            Gmpfr::get_default_precision() ) )

// CGAL_GMPFR_OP_GMPFR defines an arithmetic operation between two
// Gmpfr objects.
#define CGAL_GMPFR_OP_GMPFR(_name,_fun) \
        inline \
        Gmpfr Gmpfr::_name (const Gmpfr &a, \
                            const Gmpfr &b, \
                            std::float_round_style r){ \
                Gmpfr result(0,CGAL_GMPFR_PREC_2(a,b)); \
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

// CGAL_GMPFR_NONCOMMUTATIVE_OP defines an arithmetic operation
// between a Gmpfr object and an object of another type. All operations are
// treated as non-commutative, since MPFR does not provide all the functions
// needed to operate with any type as first operand.
#define CGAL_GMPFR_NONCOMMUTATIVE_OP(_name,_type,_member,_fun) \
        inline \
        Gmpfr Gmpfr::_name (const Gmpfr &a,_type b,std::float_round_style r){ \
                Gmpfr result(0,CGAL_GMPFR_PREC(a)); \
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

CGAL_GMPFR_OP_GMPFR(add,mpfr_add)
CGAL_GMPFR_OP_GMPFR(sub,mpfr_sub)
CGAL_GMPFR_OP_GMPFR(mul,mpfr_mul)
CGAL_GMPFR_OP_GMPFR(div,mpfr_div)

CGAL_GMPFR_NONCOMMUTATIVE_OP(add,long,b,mpfr_add_si)
CGAL_GMPFR_NONCOMMUTATIVE_OP(sub,long,b,mpfr_sub_si)
CGAL_GMPFR_NONCOMMUTATIVE_OP(mul,long,b,mpfr_mul_si)
CGAL_GMPFR_NONCOMMUTATIVE_OP(div,long,b,mpfr_div_si)

CGAL_GMPFR_NONCOMMUTATIVE_OP(add,unsigned long,b,mpfr_add_ui)
CGAL_GMPFR_NONCOMMUTATIVE_OP(sub,unsigned long,b,mpfr_sub_ui)
CGAL_GMPFR_NONCOMMUTATIVE_OP(mul,unsigned long,b,mpfr_mul_ui)
CGAL_GMPFR_NONCOMMUTATIVE_OP(div,unsigned long,b,mpfr_div_ui)

CGAL_GMPFR_NONCOMMUTATIVE_OP(add,int,b,mpfr_add_si)
CGAL_GMPFR_NONCOMMUTATIVE_OP(sub,int,b,mpfr_sub_si)
CGAL_GMPFR_NONCOMMUTATIVE_OP(mul,int,b,mpfr_mul_si)
CGAL_GMPFR_NONCOMMUTATIVE_OP(div,int,b,mpfr_div_si)

CGAL_GMPFR_NONCOMMUTATIVE_OP(add,const Gmpz&,b.mpz(),mpfr_add_z)
CGAL_GMPFR_NONCOMMUTATIVE_OP(sub,const Gmpz&,b.mpz(),mpfr_sub_z)
CGAL_GMPFR_NONCOMMUTATIVE_OP(mul,const Gmpz&,b.mpz(),mpfr_mul_z)
CGAL_GMPFR_NONCOMMUTATIVE_OP(div,const Gmpz&,b.mpz(),mpfr_div_z)

#undef CGAL_GMPFR_PREC
#undef CGAL_GMPFR_PREC_2
#undef CGAL_GMPFR_OP_GMPFR
#undef CGAL_GMPFR_NONCOMMUTATIVE_OP

// Copyright (c) 2006-2010 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_POLYNOMIAL_1_H
#define CGAL_RS_POLYNOMIAL_1_H

#include <iostream>
#include <vector>
#include <CGAL/RS/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/RS/dyadic.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <boost/operators.hpp>
#include <boost/shared_ptr.hpp>

namespace CGAL{

class Algebraic_1;
class RS_polynomial_1;
typedef boost::shared_ptr<RS_polynomial_1>        polyptr;
typedef std::pair<RS_polynomial_1,int>            polypow;
typedef std::vector<polypow>                      sqfrvec;
typedef boost::shared_ptr<sqfrvec>                sqfrptr;

class RS_polynomial_1:
        boost::addable1<RS_polynomial_1,
        boost::subtractable1<RS_polynomial_1
        > >
{
        private:
                int _capacity;
                mutable int _degree;
                mpz_t* _coef;
                mutable bool _is_sf;
                mutable polyptr _sfpart;
                mutable sqfrptr _sqfr;
                void create_storage(int);
                void free_storage();
                // fetch_gmp_functions gathers the memory functions used by
                // gmp at the object creation and stores them in _allocf,
                // _reallocf and _freef
                void *(*_allocf)(size_t);
                void *(*_reallocf)(void*,size_t,size_t);
                void (*_freef)(void*,size_t);
                void fetch_gmp_functions();
        public:
                // copy constructor and copy assignement operator
                RS_polynomial_1(const RS_polynomial_1&);
                RS_polynomial_1& operator=(const RS_polynomial_1&);
                // other constructors and destructor
                RS_polynomial_1();
                RS_polynomial_1(unsigned int);
                RS_polynomial_1(int);
                RS_polynomial_1(std::string&);
                RS_polynomial_1(mpq_srcptr);
                RS_polynomial_1(mpz_t**,int);
                ~RS_polynomial_1();
                // member functions
                void set_degree(int);
                void force_degree(int);
                int resize(int);
                void set_coef(int,mpz_srcptr);
                void set_coef(int,const CGAL::Gmpz&);
                void set_coef_ui(int,unsigned long);
                void set_coef_si(int,long);
                int get_degree()const;
                int get_degree_static()const;
                bool has_sfpart()const;
                const RS_polynomial_1& sfpart()const;
                void set_sfpart(RS_polynomial_1*)const;
                void set_sfpart(const polyptr&)const;
                void set_sf()const;
                bool has_sqfr()const;
                sqfrvec& sqfr()const;
                void set_sqfr(sqfrvec*)const;
                void set_sqfr(const sqfrptr&)const;
                mpz_ptr leading_coefficient()const;
                int first_non_zero()const;
                mpz_t* get_coefs()const;
                mpz_ptr coef(int)const;
                RS_polynomial_1& derive()const;
                RS_polynomial_1& minusx()const;
                void get_lower_bound(mpfr_ptr)const;
                void get_upper_bound(mpfr_ptr)const;
                RS_polynomial_1 times_monomial(mpz_srcptr,int)const;
                // member evaluation and sign functions
                void eval_dyadic(CGALRS_dyadic_ptr,CGALRS_dyadic_srcptr)const;
                void eval_mpfr(mpfr_ptr,mpfr_srcptr)const;
                void inexact_eval_mpfr(mpfr_ptr,mpfr_srcptr)const;
                void eval_mpfi(mpfi_ptr,mpfi_srcptr)const;
                Sign sign_dyadic(CGALRS_dyadic_srcptr)const;
                Sign sign_mpfr(mpfr_srcptr)const;
                RS::rs_sign sign_mpfi(mpfi_srcptr)const;
                double operator()(double)const;
                CGAL::Gmpz operator()(int)const;
                RS_polynomial_1 operator-()const;
                RS_polynomial_1& operator+=(const RS_polynomial_1&);
                RS_polynomial_1& operator-=(const RS_polynomial_1&);
                RS_polynomial_1 operator*(const RS_polynomial_1&)const;
                RS_polynomial_1& operator*=(const RS_polynomial_1&);
                RS_polynomial_1& operator*=(mpz_srcptr);
                RS_polynomial_1& operator*=(const CGAL::Gmpz &);
                // division is always assumed to be exact in this class
                RS_polynomial_1& operator/=(mpz_srcptr);
                RS_polynomial_1& operator/=(const CGAL::Gmpz&);
                bool operator==(const RS_polynomial_1&)const;
                // template members, including constructor
                template<class InputIt>RS_polynomial_1(InputIt,InputIt);
                template<class T>T operator()(const T&)const;
                template<class T>Sign sign_at(const T&)const;
                template<class T>RS_polynomial_1 operator*(const T&)const;
                template<class T>RS_polynomial_1& operator*=(const T&);
                template<class T>RS_polynomial_1& operator/=(const T&);
                template<class T>RS_polynomial_1 operator/(const T&)const;
};

} // namespace CGAL

#include <CGAL/RS/polynomial_1_constructors.h>
#include <CGAL/RS/polynomial_1_eval.h>
#include <CGAL/RS/polynomial_1_member.h>
#include <CGAL/RS/polynomial_1_operators.h>
#include <CGAL/RS/polynomial_1_impl.h>
#include <CGAL/RS/polynomial_1_io.h>

#endif  // CGAL_GRBS_POLYNOMIAL_1_H

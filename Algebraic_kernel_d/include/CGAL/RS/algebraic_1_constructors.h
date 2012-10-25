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

#ifndef CGAL_RS_ALGEBRAIC_1_CONSTRUCTORS_H
#define CGAL_RS_ALGEBRAIC_1_CONSTRUCTORS_H

#include <CGAL/RS/polynomial_1.h>
#include <CGAL/RS/polynomial_1_utils.h>
#include <CGAL/RS/sign_1.h>

namespace CGAL{

inline
Algebraic_1::Algebraic_1(const Algebraic_1 &i,mpfr_prec_t pl,mpfr_prec_t pr){
        if(pl>pr)
                mpfi_init2(mpfi(),pl);
        else
                mpfi_init2(mpfi(),pr);
        set_mpfi(i.mpfi());
        set_pol(i.pol());
        set_nr(i.nr());
        set_mult(i.mult());
        set_lefteval(i.lefteval());
}

inline
Algebraic_1::Algebraic_1(){
        mpfi_init(mpfi());
}

inline
Algebraic_1::Algebraic_1(int i){
        mpq_t temp;
        mpfi_init(mpfi());
        mpfi_set_si(mpfi(),(long int)i);
        mpq_init(temp);
        mpq_set_si(temp,(long int)i,1);
        RS_polynomial_1 *p=new RS_polynomial_1(temp);
        mpq_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(unsigned i){
        mpq_t temp;
        mpfi_init(mpfi());
        mpfi_set_ui(mpfi(),i);
        mpq_init(temp);
        mpq_set_ui(temp,(unsigned)i,1);
        RS_polynomial_1 *p=new RS_polynomial_1(temp);
        mpq_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(long i){
        mpq_t temp;
        mpfi_init(mpfi());
        mpfi_set_si(mpfi(),i);
        mpq_init(temp);
        mpq_set_si(temp,i,1);
        RS_polynomial_1 *p=new RS_polynomial_1(temp);
        mpq_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(unsigned long i){
        mpq_t temp;
        mpfi_init(mpfi());
        mpfi_set_ui(mpfi(),i);
        mpq_init(temp);
        mpq_set_ui(temp,i,1);
        RS_polynomial_1 *p=new RS_polynomial_1(temp);
        mpq_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(double d){
        mpq_t temp;
        mpfi_init(mpfi());
        mpfi_set_d(mpfi(),d);
        mpq_init(temp);
        mpq_set_d(temp,d);
        RS_polynomial_1 *p=new RS_polynomial_1(temp);
        mpq_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(mpz_srcptr z){
        mpq_t temp;
        mpfi_init(mpfi());
        mpfi_set_z(mpfi(),z);
        mpq_init(temp);
        mpq_set_z(temp,z);
        RS_polynomial_1 *p=new RS_polynomial_1(temp);
        mpq_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(mpq_srcptr q){
        mpfi_init(mpfi());
        mpfi_set_q(mpfi(),q);
        RS_polynomial_1 *p=new RS_polynomial_1(q);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(mpfr_srcptr src){
        Gmpfr r(src);
        mpfi_init2(mpfi(),r.get_precision());
        mpfi_set_fr(mpfi(),r.fr());
        CGAL_assertion(mpfr_equal_p(r.fr(),&mpfi()->left)!=0);
        CGAL_assertion(mpfr_equal_p(r.fr(),&mpfi()->right)!=0);
        RS_polynomial_1 *rsp=new RS_polynomial_1(Gmpq(r).mpq());
        set_pol(*rsp);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(mpfi_srcptr i){
        mpfi_init(mpfi());
        set_mpfi(i);
}

inline
Algebraic_1::Algebraic_1(const Gmpz &z){
        mpq_t temp;
        mpfi_init(mpfi());
        mpfi_set_z(mpfi(),z.mpz());
        mpq_init(temp);
        mpq_set_z(temp,z.mpz());
        RS_polynomial_1 *p=new RS_polynomial_1(temp);
        mpq_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(const Gmpq &q){
        mpfi_init(mpfi());
        mpfi_set_q(mpfi(),q.mpq());
        RS_polynomial_1 *p=new RS_polynomial_1(q.mpq());
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
Algebraic_1::Algebraic_1(const Gmpfr &r){
        mpfi_init2(mpfi(),r.get_precision());
        mpfi_set_fr(mpfi(),r.fr());
        CGAL_assertion(mpfr_equal_p(r.fr(),&mpfi()->left)!=0);
        CGAL_assertion(mpfr_equal_p(r.fr(),&mpfi()->right)!=0);
        RS_polynomial_1 *rsp=new RS_polynomial_1(Gmpq(r).mpq());
        set_pol(*rsp);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

// interesting constructor
inline
Algebraic_1::Algebraic_1(
                         mpfi_srcptr i,
                         const RS_polynomial_1 &p,
                         int n,
                         int m){
        mpfi_init(mpfi());
        set_mpfi_ptr(i);
        set_pol(p);
        set_nr(n);
        set_mult(m);
        // we don't evaluate in the sf part of p, since p is sf
        // TODO: add assertion
        set_lefteval(RSSign::signat(p,&i->left));
}

// another interesting constructor, where we have already calculated the
// left evaluation
inline
Algebraic_1::Algebraic_1(mpfi_srcptr i,const RS_polynomial_1 &p,int n,int m,
                         CGAL::Sign s){
        mpfi_init(mpfi());
        set_mpfi_ptr(i);
        set_pol(p);
        set_nr(n);
        set_mult(m);
        set_lefteval(s);
}

} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_1_CONSTRUCTORS_H

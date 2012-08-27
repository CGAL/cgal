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

#ifndef CGAL_RS_POLYNOMIAL_1_OPERATORS_H
#define CGAL_RS_POLYNOMIAL_1_OPERATORS_H

namespace CGAL{

inline
double RS_polynomial_1::operator()(double x)const{
        double result=mpz_get_d(leading_coefficient());
        int d=get_degree_static();
        for(int i=1;i<d+1;++i)
                result=x*result+mpz_get_d(coef(d-i));
        return result;
}

inline
Gmpz RS_polynomial_1::operator()(int x)const{
        Gmpz result(mpz_get_d(leading_coefficient()));
        int d=get_degree_static();
        for(int i=1;i<d+1;++i)
                result=x*result+Gmpz(coef(d-i));
        return result;
}

inline
RS_polynomial_1 RS_polynomial_1::operator-()const{
        RS_polynomial_1 opposite(_degree);
        int d=get_degree()+1;
        for(int i=0;i<d;++i)
                mpz_neg(opposite.coef(i),coef(i));
        return opposite;
}

inline
RS_polynomial_1& RS_polynomial_1::operator+=(const RS_polynomial_1 &s){
        _is_sf=false;
        _sfpart.reset();
        _sqfr.reset();
        int sd;
        if(_degree<(sd=s.get_degree())){
                mpz_t *coef_sum=(mpz_t*)(*_allocf)(sizeof(mpz_t)*(sd+1));
                for(int i=0;i<_degree+1;++i){
                        mpz_init(coef_sum[i]);
                        mpz_add(coef_sum[i],s.coef(i),coef(i));
                        mpz_clear(coef(i));
                }
                for(int i=_degree+1;i<sd+1;++i)
                        mpz_init_set(coef_sum[i],s.coef(i));
                (*_freef)(_coef,sizeof(mpz_t)*_capacity);
                _coef=coef_sum;
                _degree=sd;
        }else
                for(int i=0;i<sd+1;++i)
                        mpz_add(coef(i),s.coef(i),coef(i));
        return *this;
}

inline
RS_polynomial_1& RS_polynomial_1::operator-=(const RS_polynomial_1 &s){
        _is_sf=false;
        _sfpart.reset();
        _sqfr.reset();
        int sd;
        if(_degree<(sd=s.get_degree())){
                mpz_t *coef_sum=(mpz_t*)(*_allocf)(sizeof(mpz_t)*(sd+1));
                for(int i=0;i<_degree+1;++i){
                        mpz_init(coef_sum[i]);
                        mpz_sub(coef_sum[i],coef(i),s.coef(i));
                        mpz_clear(coef(i));
                }
                for(int i=_degree+1;i<sd+1;++i){
                        mpz_init(coef_sum[i]);
                        mpz_neg(coef_sum[i],s.coef(i));
                }
                (*_freef)(_coef,sizeof(mpz_t)*_capacity);
                _coef=coef_sum;
                _degree=sd;
        }else
                for(int i=0;i<sd+1;++i)
                        mpz_sub(coef(i),coef(i),s.coef(i));
        return *this;
}

inline
RS_polynomial_1 RS_polynomial_1::operator*(const RS_polynomial_1 &f)const{
        int d=get_degree();
        int degree_f=f.get_degree();
        int degree_p=d+degree_f;
        RS_polynomial_1 product(degree_p);
        for(int c=0;c<degree_p+1;++c){
                int max=(c<d?c:d)+1;
                for(int i=0;i<max;++i)
                        if(c-i<=degree_f)
                                mpz_addmul(product.coef(c),coef(i),f.coef(c-i));
        }
        return product;
}

inline
RS_polynomial_1& RS_polynomial_1::operator*=(const RS_polynomial_1 &f){
        _is_sf=false;
        _sfpart.reset();
        _sqfr.reset();
        RS_polynomial_1 aux(*this);
        *this=aux*f;
        return *this;
}

inline
RS_polynomial_1& RS_polynomial_1::operator*=(mpz_srcptr s){
        for(int i=0;i<_degree+1;++i)
                mpz_mul(coef(i),coef(i),s);
        return *this;
}

inline
RS_polynomial_1& RS_polynomial_1::operator*=(const CGAL::Gmpz &s){
        for(int i=0;i<_degree+1;++i)
                mpz_mul(coef(i),coef(i),s.mpz());
        return *this;
}

inline
RS_polynomial_1& RS_polynomial_1::operator/=(mpz_srcptr d){
        for(int i=0;i<_degree+1;++i)
                mpz_divexact(coef(i),coef(i),d);
        return *this;
}

inline
RS_polynomial_1& RS_polynomial_1::operator/=(const CGAL::Gmpz &d){
        for(int i=0;i<_degree+1;++i)
                mpz_divexact(coef(i),coef(i),d.mpz());
        return *this;
}

inline
bool RS_polynomial_1::operator==(const RS_polynomial_1 &p)const{
        int d=get_degree();
        int degree_p=p.get_degree();
        if(degree_p<d){
                for(int i=degree_p+1;i<d+1;++i)
                        if(mpz_sgn(coef(i)))
                                return false;
                for(int i=0;i<degree_p+1;++i)
                        if(mpz_cmp(coef(i),p.coef(i)))
                                return false;
        }else{
                for (int i=d+1;i<degree_p+1;++i)
                        if(mpz_sgn(p.coef(i)))
                                return false;
                for(int i=0;i<degree_p+1;++i)
                        if(mpz_cmp(coef(i),p.coef(i)))
                                return false;
        }
        return true;
}

} // namespace CGAL

#endif  // CGAL_RS_POLYNOMIAL_1_OPERATORS_H

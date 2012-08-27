// Copyright (c) 2009 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_POLYNOMIAL_CONVERTER
#define CGAL_RS_POLYNOMIAL_CONVERTER

#include <CGAL/Polynomial.h>
#include <CGAL/RS/polynomial_1.h>

namespace CGAL{

template<class P>
struct to_rs_poly:
public std::unary_function<P,RS_polynomial_1>{
        RS_polynomial_1 operator()(const P &p)const{
                std::cerr<<"can't convert to integer polynomial"<<std::endl;
                exit(-1);
        }
};

template<>
struct to_rs_poly<RS_polynomial_1>:
public std::unary_function<RS_polynomial_1,RS_polynomial_1>{
        const RS_polynomial_1& operator()(const RS_polynomial_1 &p)const{
                return p;
        }
};

// The conversions using this macro are efficient, since this construct an
// RS polynomial only by copying pointers. The only implementation detail
// is that it should not free the occuped by the pointed coefficients.
#define CGALRS_POLYNOMIAL_CONVERTER_REF(_T,_CONVERT) \
        template<> \
        struct to_rs_poly<Polynomial<_T> >: \
        public std::unary_function<Polynomial<_T>,RS_polynomial_1>{ \
                RS_polynomial_1& operator()(const Polynomial<_T> &p)const{ \
                        void *(*af)(size_t); \
                        void *(*rf)(void*,size_t,size_t); \
                        void (*ff)(void*,size_t); \
                        int d=p.degree(); \
                        mpz_t* c=(mpz_t*)malloc((d+1)*sizeof(mpz_t)); \
                        for(int i=0;i<=d;++i) \
                                _CONVERT; \
                        mp_get_memory_functions(&af,&rf,&ff); \
                        mp_set_memory_functions(af,rf,__cgalrs_dummy_free); \
                        RS_polynomial_1 *r=new RS_polynomial_1(&c,d); \
                        mp_set_memory_functions(af,rf,ff); \
                        return *r; \
                } \
        }

// The conversions using this macro are not intended to be efficient, since
// there is no direct way to convert from these types to mpz_t and we need
// thus to create new mpz_t's.
#define CGALRS_POLYNOMIAL_CONVERTER_COPY(_T,_CONVERT) \
        template<> \
        struct to_rs_poly<Polynomial<_T> >: \
        public std::unary_function<Polynomial<_T>,RS_polynomial_1>{ \
                RS_polynomial_1& operator()(const Polynomial<_T> &p)const{ \
                        int d=p.degree(); \
                        mpz_t* c=(mpz_t*)malloc((d+1)*sizeof(mpz_t)); \
                        for(int i=0;i<=d;++i){ \
                                mpz_init(c[i]); \
                                _CONVERT; \
                        } \
                        return *(new RS_polynomial_1(&c,d)); \
                } \
        }

//CGALRS_POLYNOMIAL_CONVERTER_REF(Gmpz,c[i][0]=*(p[i].mpz()));
CGALRS_POLYNOMIAL_CONVERTER_COPY(Gmpz,mpz_set(c[i],p[i].mpz()));
CGALRS_POLYNOMIAL_CONVERTER_COPY(int,mpz_set_si(c[i],(long)p[i]));
CGALRS_POLYNOMIAL_CONVERTER_COPY(long,mpz_set_si(c[i],p[i]));
CGALRS_POLYNOMIAL_CONVERTER_COPY(unsigned,mpz_set_ui(c[i],(unsigned long)p[i]));
CGALRS_POLYNOMIAL_CONVERTER_COPY(unsigned long,mpz_set_ui(c[i],p[i]));

#undef CGALRS_POLYNOMIAL_CONVERTER_REF
#undef CGALRS_POLYNOMIAL_CONVERTER_COPY

// convert a Gmpz rational polynomial to an integer one
template<>
struct to_rs_poly<Polynomial<Gmpq> >:
public std::unary_function<Polynomial<Gmpq>,RS_polynomial_1>{
        RS_polynomial_1& operator()(const Polynomial<Gmpq> &p)const{
                int d=p.degree();
                mpz_t denominator;
                mpz_init(denominator);
                mpz_lcm(denominator,
                        mpq_denref(p[0].mpq()),
                        mpq_denref(p[d].mpq()));
                for(int j=1;j<d;++j)
                        mpz_lcm(denominator,
                                denominator,
                                mpq_denref(p[j].mpq()));
                mpz_t* c=(mpz_t*)malloc((d+1)*sizeof(mpz_t));
                for(int i=0;i<=d;++i){
                        mpz_init(c[i]);
                        mpz_div(c[i],denominator,mpq_denref(p[i].mpq()));
                        mpz_mul(c[i],c[i],mpq_numref(p[i].mpq()));
                }
                mpz_clear(denominator);
                return *(new RS_polynomial_1(&c,d));
        }
};

template<class P>
struct from_rs_poly:
public std::unary_function<RS_polynomial_1,P>{
        P operator()(const RS_polynomial_1 &p)const{
                std::cerr<<"can't convert to integer polynomial"<<std::endl;
                exit(-1);
        }
};

template<>
struct from_rs_poly<RS_polynomial_1>:
public std::unary_function<RS_polynomial_1,RS_polynomial_1>{
        const RS_polynomial_1& operator()(const RS_polynomial_1 &p)const{
                return p;
        }
};

template<>
struct from_rs_poly<Polynomial<Gmpz> >:
public std::unary_function<RS_polynomial_1,Polynomial<Gmpz> >{
        Polynomial<Gmpz> operator()(const RS_polynomial_1 &p)const{
                typedef Polynomial_traits_d<Polynomial<Gmpz> >  PT;
                mpz_t* pcoef=p.get_coefs();
                std::vector<Gmpz> coeffs;
                int d=p.get_degree();
                for(int i=0;i<=d;++i){
                        // Gmpz c(1);
                        // *c.mpz()=*pcoef[i];
			Gmpz c(pcoef[i]);
                        coeffs.push_back(c);
                }
                return PT::Construct_polynomial()(coeffs.begin(),coeffs.end());
        }
};

template<>
struct from_rs_poly<Polynomial<Gmpq> >:
public std::unary_function<RS_polynomial_1,Polynomial<Gmpq> >{
        Polynomial<Gmpq> operator()(const RS_polynomial_1 &p)const{
                typedef Polynomial_traits_d<Polynomial<Gmpq> >  PT;
                mpz_t* pcoef=p.get_coefs();
                std::vector<Gmpq> coeffs;
                int d=p.get_degree();
                for(int i=0;i<=d;++i){
                        // Gmpq c(1);
                        // *mpq_numref(c.mpq())=*pcoef[i];
                        Gmpq c(pcoef[i]);
                        coeffs.push_back(c);
                }
                return PT::Construct_polynomial()(coeffs.begin(),coeffs.end());
        }
};

} // namespace CGAL

#endif  // CGAL_RS_POLYNOMIAL_CONVERTER

// Copyright (c) 2007-2008 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_POLYNOMIAL_FUNCTORS
#define CGAL_RS_POLYNOMIAL_FUNCTORS

#include <CGAL/basic.h>
#include <CGAL/RS/polynomial_1.h>
#include <CGAL/RS/polynomial_1_utils.h>
#include <CGAL/Exponent_vector.h>
#include <algorithm>
#include <CGAL/assertions.h>

namespace CGAL{

struct RSPolynomialFunctors{

        typedef Gmpz                    Coefficient;
        typedef Gmpz                    Innermost_coefficient;
        typedef RS_polynomial_1         Polynomial_1;
        typedef RS_polynomial_1         Type;
        static const int d=1;   // TODO: do this in a clean way...

struct Construct_polynomial{
        Polynomial_1 operator()(){
                Polynomial_1 p(2);
                p.set_coef_ui(0,0);
                p.force_degree(0);
                return p;
        };

        Polynomial_1 operator()(int i){
                Polynomial_1 p(2);
                p.set_coef_si(0,i);
                p.force_degree(0);
                return p;
        };

        template <class InputIterator>
        inline Polynomial_1 operator()(
                        InputIterator begin,InputIterator end)const{
                return Polynomial_1(begin,end);
        };

//--------------------------------------------------
//         template <class InputIterator1,class InputIterator2>
//         Polynomial_1 operator()
//         (InputIterator1 fst_coef,InputIterator1 lst_coef,
//          InputIterator2 deg)const{
//                 // the degree of the polynomial will be
//                 // the greatest degree
//                 unsigned greater=0;
//                 InputIterator1 c;
//                 InputIterator2 d=deg;
//                 for(c=fst_coef;c!=lst_coef;++c)
//                         if(d>greater)
//                                 greater=d++;
//                 // now, construct the polynomial of degree d
//                 Polynomial_1 p(d);
//                 for(c=fst_coef;c!=lst_coef;++c)
//                         p.set_coef(deg++,*c);
//                 return p;
//         };
//--------------------------------------------------

        template <class InputIterator>
        Polynomial_1 operator()
        (InputIterator begin,InputIterator end,bool is_sorted)const{
                // the degree of the polynomial will be
                // the greatest degree
                int degree;
                if(is_sorted)
                        degree=*((end-1)->first.begin());
                else{
                        degree=0;
                        for(InputIterator i=begin;i!=end;++i)
                                if(*(i->first.begin())>degree)
                                        degree=*(i->first.begin());
                }
                Polynomial_1 p(degree);
                for(InputIterator i=begin;i!=end;++i)
                        p.set_coef(*(i->first.begin()),i->second);
                return p;
        };
};      // Construct_polynomial

struct Get_coefficient:
public std::binary_function<Polynomial_1,int,Coefficient>{
        inline Coefficient operator()
                (const Polynomial_1 &p,int i)const{
                        return Coefficient(p.coef(i));
                };
};

struct Get_innermost_coefficient:
public std::binary_function<Polynomial_1,Exponent_vector,Coefficient>{
        inline Coefficient operator()(
                        const Polynomial_1 &p,Exponent_vector v)const{
                CGAL_precondition(v.size()==1);
                return Coefficient(p.coef(v[0]));
        };
};

// this is an in-place swap
struct Swap{
        inline Polynomial_1 operator()(Polynomial_1 &p,int i,int j)const{
                mpz_t temp;
                CGAL_precondition(i<=d&&j<=d);
                mpz_init_set(temp,p.coef(j));
                p.set_coef(j,p.coef(i));
                p.set_coef(i,temp);
                mpz_clear(temp);
                return p;
        };
};

struct Move{
        inline Polynomial_1& operator()(Polynomial_1 &p,int i,int j)const{
                CGAL_precondition(i<=d&&j<=d);
                p.set_coef(j,p.coef(i));
                p.set_coef_ui(i,0);
                return p;
        };
};

struct Degree{
        inline int operator()(const Polynomial_1 &p)const{
                return p.get_degree();
        };
        inline int operator()(const Polynomial_1 &p,int i)const{
                CGAL_precondition(!d);
                return p.get_degree();
        };
};

struct Total_degree:
public std::unary_function<Polynomial_1,int>{
        inline int operator()(const Polynomial_1 &p)const{
                return p.get_degree();
        };
};

struct Degree_vector:
public std::unary_function<Polynomial_1,Exponent_vector>{
        inline Exponent_vector operator()(const Polynomial_1 &p){
                Exponent_vector v(1);
                v[1]=p.get_degree();
                return v;
                // why the hell this doesn't work?
                //return Exponent_vector(1,p.get_degree());
        };
};

struct Leading_coefficient{
        inline Coefficient operator()(const Polynomial_1 &p)const{
                return Coefficient(p.leading_coefficient());
        };
        inline Coefficient operator()
                (const Polynomial_1 &p,int i)const{
                        CGAL_precondition(!i);
                        return Coefficient(p.leading_coefficient());
                };
};

struct Innermost_leading_coefficient:
public std::unary_function<Polynomial_1,Innermost_coefficient>{
        inline Innermost_coefficient operator()(const Polynomial_1 &p)const{
                return Innermost_coefficient(p.leading_coefficient());
        };
};

struct Canonicalize:
public std::unary_function<Polynomial_1,Polynomial_1>{
        Polynomial_1& operator()(Polynomial_1 &f)const{
                mpz_t temp;
                mpz_init(temp);
                Cont()(temp,f);
                f/=temp;
                mpz_clear(temp);
                return f;
        };
};

struct Derive{
        inline Polynomial_1 operator()(const Polynomial_1 &p)const{
                return p.derive();
        };
        inline Polynomial_1 operator()(const Polynomial_1 &p,int i)const{
                CGAL_precondition(!i);
                return p.derive();
        };
};

struct Evaluate{
        inline Coefficient operator()(
                        const Polynomial_1 &p,
                        const Innermost_coefficient x)const{
                return p(x);
        };
        inline Coefficient operator()(
                        const Polynomial_1 &p,
                        Innermost_coefficient c,
                        int i)const{
                CGAL_precondition(!i);
                return p(c);
        };
};

typedef Null_functor    Evaluate_homogeneous;

struct Is_zero_at{
        template <class InputIterator>
        inline bool operator()(
                        const Polynomial_1 &p,
                        InputIterator begin,
                        InputIterator end)const{
                CGAL_assertion(end=begin+1);
                return(p(*begin)==0);
        };
};

typedef Null_functor    Is_zero_at_homogeneous;

struct Sign_at{
        template <class InputIterator>
        inline Sign operator()(
                        const Polynomial_1 &p,
                        InputIterator begin,
                        InputIterator end)const{
                return p.sign_at(*begin);
        };
};

typedef Null_functor    Sign_at_homogeneous;

struct Compare:
public std::binary_function<Polynomial_1,Polynomial_1,Comparison_result>{
        Comparison_result operator()(
                        const Polynomial_1 &f,
                        const Polynomial_1 &g){
                if(f.get_degree()>g.get_degree())
                        return LARGER;
                if(f.get_degree_static()<g.get_degree_static())
                        return SMALLER;
                int i,c;
                for(i=f.get_degree_static();i>=0;--i){
                        if((c=mpz_cmp(f.coef(i),g.coef(i)))>0)
                                return LARGER;
                        if(c<0)
                                return SMALLER;
                }
                return EQUAL;
        };
};      // Compare

// now, we define five functors which are needed by the algebraic kernel,
// but they are related to polynomials (I think this is the best place
// to do it)

template <class _Gcd_policy>
struct Is_square_free_1:
        public std::unary_function<Polynomial_1,bool>{
                typedef _Gcd_policy     Gcd;
                inline bool operator()(const Polynomial_1 &p){
                        return(!(Gcd()(p,p.derive()).get_degree_static()));
                };
        };      // Is_square_free_1

template <class _Gcd_policy>
struct Make_square_free_1:
        public std::unary_function<Polynomial_1,Polynomial_1>{
                typedef _Gcd_policy     Gcd;
                inline Polynomial_1 operator()(const Polynomial_1 &p){
                        return sfpart_1<Gcd>(p);
                };
        };      // Make_square_free_1

// this function is not well defined in algebraic kernel concepts; we
// don't care, we do it correctly
template <class _Gcd_policy>
struct Square_free_factorize_1{
        template <class OutputIterator>
        int operator()(const Polynomial_1 &p,OutputIterator oi){
                typedef _Gcd_policy     Gcd;
                sqfrvec factorization(sqfr_1<Gcd>(p));
                std::copy(factorization.begin(),factorization.end(),oi);
                return factorization.size();
        };
};      // Square_free_factorize_1

template <class _Gcd_policy>
struct Is_coprime_1:
        public std::binary_function<Polynomial_1,Polynomial_1,bool>{
                inline bool operator()
                        (const Polynomial_1 &p1,const Polynomial_1 &p2){
                                typedef _Gcd_policy     Gcd;
                                return(!Gcd()(p1,p2).get_degree_static());
                        };
        };      // Is_coprime_1

template <class _Gcd_policy>
struct Make_coprime_1{
        typedef bool            result_type;
        typedef Polynomial_1    P;
        bool operator()(const P &p1,const P &p2,P &g,P &q1,P &q2){
                typedef _Gcd_policy     Gcd;
                g=Gcd()(p1,p2);
                // we don't calculate q1 and q2 when g==1, shall we?
                if(!(g.get_degree_static()))
                        return true;
                q1=*Ediv_1()(p1,g);
                q2=*Ediv_1()(p2,g);
                return false;
        };
};      // Make_coprime_1

struct IntegralDivision:
        public std::binary_function<Polynomial_1,Polynomial_1,Polynomial_1>{
                inline Polynomial_1 operator()(
                                const Polynomial_1 &p1,
                                const Polynomial_1 &p2){
                        return *Ediv_1()(p1,p2);
                };
};      // IntegralDivision

typedef IntegralDivision        IntegralDivisionUpToConstantFactor;

};      // RSPolynomialFunctors

} // namespace CGAL

#endif  // CGAL_RS_POLYNOMIAL_FUNCTORS

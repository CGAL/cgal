// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.

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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_ALGEBRAIC_1_H
#define CGAL_RS_ALGEBRAIC_1_H

#include <boost/operators.hpp>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Gmpq.h>
#include <iostream>

namespace CGAL{
namespace RS_AK1{

// This class represents the simplest algebraic number one can think about.
// One algebraic number is represented by the polynomial of which it is
// root and the two endpoints of an interval that contains it, and no other
// root. Polynomial type and bound type are the first two template
// parameters.
//
// The third template parameter is a refiner, a function object that
// receives the polynomial and the bounds defining an algebraic number and
// an integer p, and modifies the two bounds until the difference between
// the two bounds is less than x*2^(-p), where x is the value of the
// represented algebraic number. The signature of a refiner must be:
//      void
//      Refiner_()(const Polynomial_&,Bound_&,Bound_&,int p);
//
// The fourth template argument is a comparator, a function object that
// receives the polynomials and bounds defining two algebraic numbres and
// just compares them, returning a CGAL::Comparison_result. The signature
// of a comparator must be:
//      CGAL::Comparison_result
//      Comparator_()(const Polynomial_&,Bound_&,Bound_&,
//                    const Polynomial_&,Bound_&,Bound_&);
// The comparator can modify the bounds, with the condition that the
// algebraic numbers remain consistent (one and only one root on each
// interval).

template <class Polynomial_,
          class Bound_,
          class Refiner_,
          class Comparator_,
          class Ptraits_>
class Algebraic_1:
boost::totally_ordered<Algebraic_1<Polynomial_,
                                   Bound_,
                                   Refiner_,
                                   Comparator_,
                                   Ptraits_>,
                       double>{
        protected:
        typedef Polynomial_                             Polynomial;
        typedef Bound_                                  Bound;
        typedef Refiner_                                Refiner;
        typedef Comparator_                             Comparator;
        typedef Ptraits_                                Ptraits;
        typedef typename Ptraits::Coefficient_type      Coefficient;
        typedef typename Ptraits::Scale                 Scale;
        typedef Algebraic_1<Polynomial,
                            Bound,
                            Refiner,
                            Comparator,
                            Ptraits>
                                                        Algebraic;

        Polynomial pol;
        mutable Bound left,right;

        public:
        Algebraic_1(){};
        Algebraic_1(const Polynomial &p,
                    const Bound &l,
                    const Bound &r):pol(p),left(l),right(r){
                CGAL_assertion(l<=r);
        }
        Algebraic_1(const Algebraic &a):
        pol(a.pol),left(a.left),right(a.right){}
        // XXX: This assumes that Gmpq is constructible from bound type and
        // that polynomial coefficient type is constructible from mpz_t.
        Algebraic_1(const Bound &b):left(b),right(b){
                typedef typename Ptraits::Shift         shift;
                Gmpq q(b);
                pol=Coefficient(mpq_denref(q.mpq()))*
                        shift()(Polynomial(1),1,0)-
                        Coefficient(mpq_numref(q.mpq()));
                CGAL_assertion(left==right&&left==b);
        }
        // XXX: This implementation assumes that the bound type is Gmpfr
        // and that T can be exactly converted to Gmpq. This constructor
        // can handle types such as int, unsigned and long. 
        template <class T>
        Algebraic_1(const T &t){
                typedef typename Ptraits::Shift         shift;
                CGAL::Gmpq q(t);
                pol=Coefficient(mpq_denref(q.mpq()))*
                    shift()(Polynomial(1),1,0)-
                        Coefficient(mpq_numref(q.mpq()));
                left=Bound(t,std::round_toward_neg_infinity);
                right=Bound(t,std::round_toward_infinity);
                CGAL_assertion(left<=t&&right>=t);
        }
        // XXX: This constructor assumes the bound is a Gmpfr.
        Algebraic_1(const CGAL::Gmpq &q){
                typedef typename Ptraits::Shift         shift;
                pol=Coefficient(mpq_denref(q.mpq()))*
                    shift()(Polynomial(1),1,0)-
                        Coefficient(mpq_numref(q.mpq()));
                left=Bound();
                right=Bound();
                mpfr_t b;
                mpfr_init(b);
                mpfr_set_q(b,q.mpq(),GMP_RNDD);
                mpfr_swap(b,left.fr());
                mpfr_set_q(b,q.mpq(),GMP_RNDU);
                mpfr_swap(b,right.fr());
                mpfr_clear(b);
                CGAL_assertion(left<=q&&right>=q);
        }
        ~Algebraic_1(){}

        Algebraic_1& operator=(const Algebraic &a){
                pol=a.get_pol();
                left=a.get_left();
                right=a.get_right();
                return *this;
        }

        Polynomial get_pol()const{return pol;}
        Bound& get_left()const{return left;}
        Bound& get_right()const{return right;}

        Algebraic operator-()const{
                return Algebraic(Scale()(get_pol(),Coefficient(-1)),
                                 -right,
                                 -left);
        }

#define CGAL_RS_COMPARE_ALGEBRAIC(_a) \
        (Comparator()(get_pol(),get_left(),get_right(), \
                      (_a).get_pol(),(_a).get_left(),(_a).get_right()))

        Comparison_result compare(Algebraic a)const{
                return CGAL_RS_COMPARE_ALGEBRAIC(a);
        };

#define CGAL_RS_COMPARE_ALGEBRAIC_TYPE(_t) \
        bool operator<(_t t)const \
        {Algebraic a(t);return CGAL_RS_COMPARE_ALGEBRAIC(a)==CGAL::SMALLER;} \
        bool operator>(_t t)const \
        {Algebraic a(t);return CGAL_RS_COMPARE_ALGEBRAIC(a)==CGAL::LARGER;} \
        bool operator==(_t t)const \
        {Algebraic a(t);return CGAL_RS_COMPARE_ALGEBRAIC(a)==CGAL::EQUAL;}

        bool operator==(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC(a)==CGAL::EQUAL;}
        bool operator!=(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC(a)!=CGAL::EQUAL;}
        bool operator<(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC(a)==CGAL::SMALLER;}
        bool operator<=(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC(a)!=CGAL::LARGER;}
        bool operator>(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC(a)==CGAL::LARGER;}
        bool operator>=(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC(a)!=CGAL::SMALLER;}

        CGAL_RS_COMPARE_ALGEBRAIC_TYPE(double)

#undef CGAL_RS_COMPARE_ALGEBRAIC_TYPE
#undef CGAL_RS_COMPARE_ALGEBRAIC

#ifdef IEEE_DBL_MANT_DIG
#define CGAL_RS_DBL_PREC        IEEE_DBL_MANT_DIG
#else
#define CGAL_RS_DBL_PREC        53
#endif
        // This function is const because left and right are mutable.
        double to_double()const{
                typedef Real_embeddable_traits<Bound>                   RT;
                typedef typename RT::To_double                          TD;
                Refiner()(pol,left,right,CGAL_RS_DBL_PREC);
                CGAL_assertion(TD()(left)==TD()(right));
                return TD()(left);
        }
        std::pair<double,double> to_interval()const{
                typedef Real_embeddable_traits<Bound>                   RT;
                typedef typename RT::To_interval                        TI;
                return std::make_pair(TI()(get_left()).first,
                                      TI()(get_right()).second);
        }
#undef CGAL_RS_DBL_PREC

        void set_left(const Bound &l)const{
                left=l;
        }
        void set_right(const Bound &r)const{
                right=r;
        }
        void set_pol(const Polynomial &p){
                pol=p;
        }

}; // class Algebraic_1

} // namespace RS_AK1

// We define Algebraic_1 as real-embeddable
template <class Polynomial_,
          class Bound_,
          class Refiner_,
          class Comparator_,
          class Ptraits_>
class Real_embeddable_traits<RS_AK1::Algebraic_1<Polynomial_,
                                                 Bound_,
                                                 Refiner_,
                                                 Comparator_,
                                                 Ptraits_> >:
public INTERN_RET::Real_embeddable_traits_base<
                RS_AK1::Algebraic_1<Polynomial_,
                                    Bound_,
                                    Refiner_,
                                    Comparator_,
                                    Ptraits_>,
                CGAL::Tag_true>{
        typedef Polynomial_                             P;
        typedef Bound_                                  B;
        typedef Refiner_                                R;
        typedef Comparator_                             C;
        typedef Ptraits_                                T;

        public:

        typedef RS_AK1::Algebraic_1<P,B,R,C,T> Type;
        typedef CGAL::Tag_true                          Is_real_embeddable;
        typedef bool                                    Boolean;
        typedef CGAL::Sign                              Sign;
        typedef CGAL::Comparison_result                 Comparison_result;

        typedef INTERN_RET::Real_embeddable_traits_base<Type,CGAL::Tag_true>
                                                        Base;
        typedef typename Base::Compare                  Compare;

        class Sgn:public CGAL::unary_function<Type,CGAL::Sign>{
                public:
                CGAL::Sign operator()(const Type &a)const{
                        return Compare()(a,Type(0));
                }
        };

        class To_double:public CGAL::unary_function<Type,double>{
                public:
                double operator()(const Type &a)const{return a.to_double();}
        };

        class To_interval:
        public CGAL::unary_function<Type,std::pair<double,double> >{
                public:
                std::pair<double,double> operator()(const Type &a)const{
                        return a.to_interval();
                }
        };

        class Is_zero:public CGAL::unary_function<Type,Boolean>{
                public:
                bool operator()(const Type &a)const{
                        return Sgn()(a)==CGAL::ZERO;
                }
        };

        class Is_finite:public CGAL::unary_function<Type,Boolean>{
                public:
                bool operator()(const Type&)const{return true;}
        };

        class Abs:public CGAL::unary_function<Type,Type>{
                public:
                Type operator()(const Type &a)const{
                        return Sgn()(a)==CGAL::NEGATIVE?-a:a;
                }
        };
};

template <class P,class B,class R,class C,class T>
inline
RS_AK1::Algebraic_1<P,B,R,C,T> min
BOOST_PREVENT_MACRO_SUBSTITUTION(RS_AK1::Algebraic_1<P,B,R,C,T> a,
                                 RS_AK1::Algebraic_1<P,B,R,C,T> b){
        return(a<b?a:b);
}

template <class P,class B,class R,class C,class T>
inline
RS_AK1::Algebraic_1<P,B,R,C,T> max
BOOST_PREVENT_MACRO_SUBSTITUTION(RS_AK1::Algebraic_1<P,B,R,C,T> a,
                                 RS_AK1::Algebraic_1<P,B,R,C,T> b){
        return(a>b?a:b);
}

template <class P,class B,class R,class C,class T>
inline
std::ostream& operator<<(std::ostream &o,
                         const RS_AK1::Algebraic_1<P,B,R,C,T> &a){
        return(o<<'['<<a.get_pol()<<','<<
               a.get_left()<<','<<
               a.get_right()<<']');
}

// XXX: This function works, but it will be nice to rewrite it cleanly.
template <class P,class B,class R,class C,class T>
inline
std::istream& operator>>(std::istream &i,
                         RS_AK1::Algebraic_1<P,B,R,C,T> &a){
        std::istream::int_type c;
        P pol;
        B lb,rb;
        c=i.get();
        if(c!='['){
                CGAL_error_msg("error reading istream, \'[\' expected");
                return i;
        }
        i>>pol;
        c=i.get();
        if(c!=','){
                CGAL_error_msg("error reading istream, \',\' expected");
                return i;
        }
        i>>lb;
        c=i.get();
        if(c!=','){
                CGAL_error_msg("error reading istream, \',\' expected");
                return i;
        }
        i>>rb;
        c=i.get();
        if(c!=']'){
                CGAL_error_msg("error reading istream, \']\' expected");
                return i;
        }
        a=RS_AK1::Algebraic_1<P,B,R,C,T>(pol,lb,rb);
        return i;
}

} // namespace CGAL

#endif // CGAL_RS_ALGEBRAIC_1_H

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
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_ALGEBRAIC_Z_1_H
#define CGAL_RS_ALGEBRAIC_Z_1_H

#include "algebraic_1.h"

namespace CGAL{
namespace RS_AK1{

// This class represents an algebraic number storing two polynomials of
// which it is root: one of them is the given polynomial; the other one is
// an integer polynomial, which is used to perform the operations such as
// refinements. Since RS works only with integer polynomials, this
// architecture permits to convert the input polynomials only once.

template <class Polynomial_,
          class ZPolynomial_,
          class Bound_,
          class ZRefiner_,
          class ZComparator_,
          class Ptraits_,
          class ZPtraits_>
class Algebraic_z_1:
public Algebraic_1<Polynomial_,Bound_,ZRefiner_,ZComparator_,ZPtraits_>,
boost::totally_ordered<Algebraic_z_1<Polynomial_,
                                     ZPolynomial_,
                                     Bound_,
                                     ZRefiner_,
                                     ZComparator_,
                                     Ptraits_,
                                     ZPtraits_>,
                       double>{
        protected:
        typedef Polynomial_                             Polynomial;
        typedef ZPolynomial_                            ZPolynomial;
        typedef Bound_                                  Bound;
        typedef ZRefiner_                               ZRefiner;
        typedef ZComparator_                            ZComparator;
        typedef ZPtraits_                               ZPtraits;
        typedef typename ZPtraits::Coefficient_type     ZCoefficient;
        typedef typename ZPtraits::Scale                ZScale;
        typedef Ptraits_                                Ptraits;
        typedef typename Ptraits::Coefficient_type      Coefficient;
        typedef typename Ptraits::Scale                 Scale;
        typedef Algebraic_1<Polynomial,Bound,ZRefiner,ZComparator,ZPtraits>
                                                        Algebraic_base;
        typedef Algebraic_z_1<Polynomial,
                              ZPolynomial,
                              Bound,
                              ZRefiner,
                              ZComparator,
                              Ptraits,
                              ZPtraits>
                                                        Algebraic;
        ZPolynomial zpol;

        public:
        Algebraic_z_1(){};
        Algebraic_z_1(const Polynomial &p,
                      const ZPolynomial &zp,
                      const Bound &l,
                      const Bound &r):Algebraic_base(p,l,r),zpol(zp){
                CGAL_assertion(l<=r);
        }
        Algebraic_z_1(const Algebraic &a):
        Algebraic_base(a.pol,a.left,a.right),zpol(a.zpol){}
        Algebraic_z_1(const Bound &b){
                typedef typename Ptraits::Shift         Shift;
                typedef typename ZPtraits::Shift        ZShift;
                this->left=b;
                this->right=b;
                Gmpq q(b);
                this->pol=Coefficient(mpq_denref(q.mpq()))*
                        Shift()(Polynomial(1),1,0)-
                        Coefficient(mpq_numref(q.mpq()));
                zpol=ZCoefficient(mpq_denref(q.mpq()))*
                        ZShift()(ZPolynomial(1),1,0)-
                        ZCoefficient(mpq_numref(q.mpq()));
                CGAL_assertion(this->left==this->right);
                CGAL_assertion(this->left==b);
        }
        template <class T>
        Algebraic_z_1(const T &t){
                typedef typename Ptraits::Shift         Shift;
                typedef typename ZPtraits::Shift        ZShift;
                CGAL::Gmpq q(t);
                this->pol=Coefficient(mpq_denref(q.mpq()))*
                    Shift()(Polynomial(1),1,0)-
                        Coefficient(mpq_numref(q.mpq()));
                zpol=ZCoefficient(mpq_denref(q.mpq()))*
                    ZShift()(ZPolynomial(1),1,0)-
                        ZCoefficient(mpq_numref(q.mpq()));
                this->left=Bound(t,std::round_toward_neg_infinity);
                this->right=Bound(t,std::round_toward_infinity);
                CGAL_assertion(this->left<=t&&this->right>=t);
        }
        Algebraic_z_1(const CGAL::Gmpq &q){
                typedef typename Ptraits::Shift         Shift;
                typedef typename ZPtraits::Shift        ZShift;
                this->pol=Coefficient(mpq_denref(q.mpq()))*
                    Shift()(Polynomial(1),1,0)-
                        Coefficient(mpq_numref(q.mpq()));
                zpol=ZCoefficient(mpq_denref(q.mpq()))*
                    ZShift()(ZPolynomial(1),1,0)-
                        ZCoefficient(mpq_numref(q.mpq()));
                this->left=Bound();
                this->right=Bound();
                mpfr_t b;
                mpfr_init(b);
                mpfr_set_q(b,q.mpq(),GMP_RNDD);
                mpfr_swap(b,this->left.fr());
                mpfr_set_q(b,q.mpq(),GMP_RNDU);
                mpfr_swap(b,this->right.fr());
                mpfr_clear(b);
                CGAL_assertion(this->left<=q&&this->right>=q);
        }
        ~Algebraic_z_1(){}

        ZPolynomial get_zpol()const{return zpol;}
        void set_zpol(const ZPolynomial &p){zpol=p;}

        Algebraic operator-()const{
                return Algebraic(Scale()(this->get_pol(),Coefficient(-1)),
                                 ZScale()(get_zpol(),ZCoefficient(-1)),
                                 -this->right,
                                 -this->left);
        }

#define CGAL_RS_COMPARE_ALGEBRAIC_Z(_a) \
        (ZComparator()(get_zpol(),this->get_left(),this->get_right(), \
                      (_a).get_zpol(),(_a).get_left(),(_a).get_right()))

        Comparison_result compare(Algebraic a)const{
                return CGAL_RS_COMPARE_ALGEBRAIC_Z(a);
        };

#define CGAL_RS_COMPARE_ALGEBRAIC_Z_TYPE(_t) \
        bool operator<(_t t)const \
        {Algebraic a(t);return CGAL_RS_COMPARE_ALGEBRAIC_Z(a)==CGAL::SMALLER;} \
        bool operator>(_t t)const \
        {Algebraic a(t);return CGAL_RS_COMPARE_ALGEBRAIC_Z(a)==CGAL::LARGER;} \
        bool operator==(_t t)const \
        {Algebraic a(t);return CGAL_RS_COMPARE_ALGEBRAIC_Z(a)==CGAL::EQUAL;}

        bool operator==(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC_Z(a)==CGAL::EQUAL;}
        bool operator!=(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC_Z(a)!=CGAL::EQUAL;}
        bool operator<(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC_Z(a)==CGAL::SMALLER;}
        bool operator<=(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC_Z(a)!=CGAL::LARGER;}
        bool operator>(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC_Z(a)==CGAL::LARGER;}
        bool operator>=(Algebraic a)const
                {return CGAL_RS_COMPARE_ALGEBRAIC_Z(a)!=CGAL::SMALLER;}

        CGAL_RS_COMPARE_ALGEBRAIC_Z_TYPE(double)

#undef CGAL_RS_COMPARE_ALGEBRAIC_Z_TYPE
#undef CGAL_RS_COMPARE_ALGEBRAIC_Z

#ifdef IEEE_DBL_MANT_DIG
#define CGAL_RS_DBL_PREC        IEEE_DBL_MANT_DIG
#else
#define CGAL_RS_DBL_PREC        53
#endif
        double to_double()const{
                typedef Real_embeddable_traits<Bound>                   RT;
                typedef typename RT::To_double                          TD;
                ZRefiner()(get_zpol(),
                           this->get_left(),
                           this->get_right(),
                           CGAL_RS_DBL_PREC);
                CGAL_assertion(TD()(this->get_left())==TD()(this->get_right()));
                return TD()(this->get_left());
        }
        std::pair<double,double> to_interval()const{
                typedef Real_embeddable_traits<Bound>                   RT;
                typedef typename RT::To_interval                        TI;
                return std::make_pair(TI()(this->get_left()).first,
                                      TI()(this->get_right()).second);
        }
#undef CGAL_RS_DBL_PREC
}; // class Algebraic_z_1

} // namespace RS_AK1

// We define Algebraic_z_1 as real-embeddable
template <class Polynomial_,
          class ZPolynomial_,
          class Bound_,
          class ZRefiner_,
          class ZComparator_,
          class Ptraits_,
          class ZPtraits_>
class Real_embeddable_traits<RS_AK1::Algebraic_z_1<Polynomial_,
                                                   ZPolynomial_,
                                                   Bound_,
                                                   ZRefiner_,
                                                   ZComparator_,
                                                   Ptraits_,
                                                   ZPtraits_> >:
public INTERN_RET::Real_embeddable_traits_base<
                RS_AK1::Algebraic_z_1<Polynomial_,
                                      ZPolynomial_,
                                      Bound_,
                                      ZRefiner_,
                                      ZComparator_,
                                      Ptraits_,
                                      ZPtraits_>,
                CGAL::Tag_true>{
        typedef Polynomial_                             P;
        typedef ZPolynomial_                            ZP;
        typedef Bound_                                  B;
        typedef ZRefiner_                               R;
        typedef ZComparator_                            C;
        typedef Ptraits_                                T;
        typedef ZPtraits_                               ZT;

        public:

        typedef RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT>  Type;
        typedef CGAL::Tag_true                          Is_real_embeddable;
        typedef bool                                    Boolean;
        typedef CGAL::Sign                              Sign;
        typedef CGAL::Comparison_result                 Comparison_result;

        typedef INTERN_RET::Real_embeddable_traits_base<Type,CGAL::Tag_true>
                                                        Base;
        typedef typename Base::Compare                  Compare;

        class Sgn:public std::unary_function<Type,CGAL::Sign>{
                public:
                CGAL::Sign operator()(const Type &a)const{
                        return Compare()(a,Type(0));
                }
        };

        class To_double:public std::unary_function<Type,double>{
                public:
                double operator()(const Type &a)const{return a.to_double();}
        };

        class To_interval:
        public std::unary_function<Type,std::pair<double,double> >{
                public:
                std::pair<double,double> operator()(const Type &a)const{
                        return a.to_interval();
                }
        };

        class Is_zero:public std::unary_function<Type,Boolean>{
                public:
                bool operator()(const Type &a)const{
                        return Sgn()(a)==CGAL::ZERO;
                }
        };

        class Is_finite:public std::unary_function<Type,Boolean>{
                public:
                bool operator()(const Type&)const{return true;}
        };

        class Abs:public std::unary_function<Type,Type>{
                public:
                Type operator()(const Type &a)const{
                        return Sgn()(a)==CGAL::NEGATIVE?-a:a;
                }
        };
};

template <class P,class ZP,class B,class R,class C,class T,class ZT>
inline
RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT> min
BOOST_PREVENT_MACRO_SUBSTITUTION(RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT> a,
                                 RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT> b){
        return(a<b?a:b);
}

template <class P,class ZP,class B,class R,class C,class T,class ZT>
inline
RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT> max
BOOST_PREVENT_MACRO_SUBSTITUTION(RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT> a,
                                 RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT> b){
        return(a>b?a:b);
}

template <class P,class ZP,class B,class R,class C,class T,class ZT>
inline
std::ostream& operator<<(std::ostream &o,
                         const RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT> &a){
        return(o<<'['<<a.get_pol()<<','<<
               a.get_zpol()<<','<<
               a.get_left()<<','<<
               a.get_right()<<']');
}

template <class P,class ZP,class B,class R,class C,class T,class ZT>
inline
std::istream& operator>>(std::istream &i,
                         RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT> &a){
        std::istream::int_type c;
        P pol;
        ZP zpol;
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
        i>>zpol;
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
        a=RS_AK1::Algebraic_z_1<P,ZP,B,R,C,T,ZT>(pol,zpol,lb,rb);
        return i;
}

} // namespace CGAL

#endif // CGAL_RS_ALGEBRAIC_Z_1_H

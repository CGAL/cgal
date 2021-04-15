// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_FUNCTORS_Z_1_H
#define CGAL_RS_FUNCTORS_Z_1_H

#include <vector>
#include <CGAL/Gmpfi.h>

namespace CGAL{
namespace RS_AK1{

template <class Polynomial_,
          class ZPolynomial_,
          class PolConverter_,
          class Algebraic_,
          class Bound_,
          class Coefficient_,
          class Isolator_>
struct Construct_algebraic_real_z_1{
        typedef Polynomial_                                     Polynomial;
        typedef ZPolynomial_                                    ZPolynomial;
        typedef PolConverter_                                   PolConverter;
        typedef Algebraic_                                      Algebraic;
        typedef Bound_                                          Bound;
        typedef Coefficient_                                    Coefficient;
        typedef Isolator_                                       Isolator;
        typedef Algebraic                                       result_type;

        template <class T>
        Algebraic operator()(const T &a)const{
                return Algebraic(a);
        }

        Algebraic operator()(const Polynomial &p,size_t i)const{
                ZPolynomial zp=PolConverter()(p);
                Isolator isol(zp);
                return Algebraic(p,
                                 zp,
                                 isol.left_bound(i),
                                 isol.right_bound(i));
        }

        Algebraic operator()(const Polynomial &p,
                             const Bound &l,
                             const Bound &r)const{
                return Algebraic(p,PolConverter()(p),l,r);
        }

}; // struct Construct_algebraic_real_z_1

template <class Polynomial_,class Algebraic_>
struct Compute_polynomial_z_1:
public CGAL::cpp98::unary_function<Algebraic_,Polynomial_>{
        typedef Polynomial_                                     Polynomial;
        typedef Algebraic_                                      Algebraic;
        Polynomial operator()(const Algebraic &x)const{
                return x.get_pol();
        }
}; // struct Compute_polynomial_z_1

template <class Polynomial_,class Ptraits_>
struct Is_coprime_z_1:
public CGAL::cpp98::binary_function<Polynomial_,Polynomial_,bool>{
        typedef Polynomial_                                     Polynomial;
        typedef Ptraits_                                        Ptraits;
        typedef typename Ptraits::Gcd_up_to_constant_factor     Gcd;
        typedef typename Ptraits::Degree                        Degree;
        inline bool operator()(const Polynomial &p1,const Polynomial &p2)const{
                return Degree()(Gcd()(p1,p2))==0;
        }
}; // struct Is_coprime_z_1

template <class Polynomial_,class Ptraits_>
struct Make_coprime_z_1{
        typedef Polynomial_                                     Polynomial;
        typedef Ptraits_                                        Ptraits;
        typedef typename Ptraits::Gcd_up_to_constant_factor     Gcd;
        typedef typename Ptraits::Degree                        Degree;
        typedef typename Ptraits::Integral_division_up_to_constant_factor
                                                                IDiv;
        bool operator()(const Polynomial &p1,
                        const Polynomial &p2,
                        Polynomial &g,
                        Polynomial &q1,
                        Polynomial &q2)const{
                g=Gcd()(p1,p2);
                q1=IDiv()(p1,g);
                q2=IDiv()(p2,g);
                return Degree()(Gcd()(p1,p2))==0;
        }
}; // struct Make_coprime_z_1

template <class Polynomial_,
          class ZPolynomial_,
          class PolConverter_,
          class Bound_,
          class Algebraic_,
          class Isolator_,
          class Signat_,
          class Ptraits_,
          class ZPtraits_>
struct Solve_z_1{
        typedef Polynomial_                                     Polynomial_1;
        typedef ZPolynomial_                                    ZPolynomial_1;
        typedef PolConverter_                                   PolConverter;
        typedef Bound_                                          Bound;
        typedef Algebraic_                                      Algebraic;
        typedef Isolator_                                       Isolator;
        typedef Signat_                                         ZSignat;
        typedef Ptraits_                                        Ptraits;
        typedef typename Ptraits::Gcd_up_to_constant_factor     Gcd;
        typedef typename Ptraits::Square_free_factorize_up_to_constant_factor
                                                                Sqfr;
        typedef typename Ptraits::Degree                        Degree;
        typedef typename Ptraits::Make_square_free              Sfpart;
        typedef ZPtraits_                                       ZPtraits;
        typedef typename ZPtraits::Gcd_up_to_constant_factor    ZGcd;
        typedef typename ZPtraits::Square_free_factorize_up_to_constant_factor
                                                                ZSqfr;
        typedef typename ZPtraits::Degree                       ZDegree;
        typedef typename ZPtraits::Make_square_free             ZSfpart;

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  OutputIterator res)const{
                typedef std::pair<ZPolynomial_1,int>    zpolmult;
                typedef std::vector<zpolmult>           zsqvec;

                ZPolynomial_1 zp=PolConverter()(p);
                Polynomial_1 sfp=Sfpart()(p);
                ZPolynomial_1 zsfp=PolConverter()(sfp);
                zsqvec zsfv;
                ZSqfr()(zp,std::back_inserter(zsfv));
                Isolator isol(zsfp);
                int *m=(int*)calloc(isol.number_of_real_roots(),sizeof(int));
                for(typename zsqvec::iterator i=zsfv.begin();i!=zsfv.end();++i){
                        int k=ZDegree()(i->first);
                        ZSignat signof(i->first);
                        for(int j=0;k&&j<isol.number_of_real_roots();++j){
                                if(!m[j]){
                                        CGAL::Sign sg_l=
                                                signof(isol.left_bound(j));
                                        CGAL::Sign sg_r=
                                                signof(isol.right_bound(j));
                                        if((sg_l!=sg_r)||
                                           ((sg_l==CGAL::ZERO)&&
                                            (sg_r==CGAL::ZERO))){
                                                m[j]=i->second;
                                                --k;
                                        }
                                }
                        }
                }
                for(int l=0;l<isol.number_of_real_roots();++l)
                        *res++=std::make_pair(Algebraic(sfp,
                                                        zsfp,
                                                        isol.left_bound(l),
                                                        isol.right_bound(l)),
                                              m[l]);
                free(m);
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  bool,
                                  OutputIterator res)const{
                ZPolynomial_1 zp=PolConverter()(p);
                Isolator isol(zp);
                for(int l=0;l<isol.number_of_real_roots();++l)
                        *res++=Algebraic(p,
                                         zp,
                                         isol.left_bound(l),
                                         isol.right_bound(l));
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  const Bound &l,
                                  const Bound &u,
                                  OutputIterator res)const{
                typedef std::vector<Algebraic>                  RV;
                typedef std::pair<Polynomial_1,int>             PM;
                typedef std::vector<PM>                         PMV;
                typedef typename PMV::iterator                  PMVI;
                CGAL_precondition_msg(l<=u,
                                      "left bound must be <= right bound");
                RV roots; // all roots of the polynomial
                this->operator()(p,false,std::back_inserter(roots));
                size_t nb_roots=roots.size();
                // indices of the first and last roots to be reported:
                size_t index_l=0,index_u;
                while(index_l<nb_roots&&roots[index_l]<l)
                        ++index_l;
                CGAL_assertion(index_l<=nb_roots);
                if(index_l==nb_roots)
                        return res;
                index_u=index_l;
                while(index_u<nb_roots&&roots[index_u]<u)
                        ++index_u;
                CGAL_assertion(index_u<=nb_roots);
                if(index_u==index_l)
                        return res;
                // now, we have to return roots in [index_l,index_u)
                PMV sfv;
                Sqfr()(p,std::back_inserter(sfv)); // square-free fact. of p
                // array to store the multiplicities
                int *m=(int*)calloc(nb_roots,sizeof(int));
                // we iterate over all the pairs <root,mult> and match the
                // roots in the interval [index_l,index_u)
                for(PMVI i=sfv.begin();i!=sfv.end();++i){
                        ZPolynomial_1 zifirst=PolConverter()(i->first);
                        int k=ZDegree()(zifirst);
                        ZSignat signof(zifirst);
                        for(size_t j=index_l;k&&j<index_u;++j){
                                if(!m[j]){
                                        CGAL::Sign sg_l=
                                                signof(roots[j].get_left());
                                        CGAL::Sign sg_r=
                                                signof(roots[j].get_right());
                                        if((sg_l!=sg_r)||
                                           ((sg_l==CGAL::ZERO)&&
                                            (sg_r==CGAL::ZERO))){
                                                m[j]=i->second;
                                                --k;
                                        }
                                }
                        }
                }
                for(size_t l=index_l;l<index_u;++l)
                        *res++=std::make_pair(roots[l],m[l]);
                free(m);
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  bool known_to_be_square_free,
                                  const Bound &l,
                                  const Bound &u,
                                  OutputIterator res)const{
                typedef std::vector<Algebraic>                  RV;
                typedef typename RV::iterator                   RVI;
                CGAL_precondition_msg(l<=u,
                                      "left bound must be <= right bound");
                RV roots;
                this->operator()(p,
                                 known_to_be_square_free,
                                 std::back_inserter(roots));
                for(RVI it=roots.begin();it!=roots.end();it++)
                        if(*it>=l&&*it<=u)
                                *res++=*it;
                return res;
        }

}; // struct Solve_z_1

template <class Polynomial_,
          class ZPolynomial_,
          class PolConverter_,
          class Bound_,
          class Algebraic_,
          class Refiner_,
          class Signat_,
          class Ptraits_,
          class ZPtraits_>
class Sign_at_z_1:
public CGAL::cpp98::binary_function<Polynomial_,Algebraic_,CGAL::Sign>{
        // This implementation will work with any polynomial type whose
        // coefficient type is explicit interoperable with Gmpfi.
        // TODO: Make this function generic.
        public:
        typedef Polynomial_                                     Polynomial_1;
        typedef ZPolynomial_                                    ZPolynomial_1;
        typedef PolConverter_                                   PolConverter;
        typedef Bound_                                          Bound;
        typedef Algebraic_                                      Algebraic;
        typedef Refiner_                                        Refiner;
        typedef Signat_                                         ZSignat;
        typedef Ptraits_                                        Ptraits;
        typedef ZPtraits_                                       ZPtraits;

        private:
        CGAL::Uncertain<CGAL::Sign> eval_interv(const Polynomial_1 &p,
                                                const Bound &l,
                                                const Bound &r)const{
                typedef typename Ptraits::Substitute                    Subst;
                std::vector<CGAL::Gmpfi> substitutions;
                substitutions.push_back(CGAL::Gmpfi(l,r));
                CGAL::Gmpfi eval=Subst()(p,
                                         substitutions.begin(),
                                         substitutions.end());
                return eval.sign();
        }

        // This function assumes that the sign of the evaluation is not zero,
        // it just refines x until having the correct sign.
        CGAL::Sign refine_and_return(const Polynomial_1 &p,Algebraic x)const{
                CGAL::Gmpfr xl(x.get_left());
                CGAL::Gmpfr xr(x.get_right());
                CGAL::Uncertain<CGAL::Sign> s;
                for(;;){
                        Refiner()(x.get_zpol(),
                                  xl,
                                  xr,
                                  2*CGAL::max(xl.get_precision(),
                                              xr.get_precision()));
                        s=eval_interv(p,xl,xr);
                        if(!s.is_same(Uncertain<CGAL::Sign>::indeterminate())){
                                x.set_left(xl);
                                x.set_right(xr);
                                return s;
                        }
                }
        }

        public:
        CGAL::Sign operator()(const Polynomial_1 &p,Algebraic x)const{
                typedef typename Ptraits::Gcd_up_to_constant_factor     Gcd;
                typedef typename Ptraits::Make_square_free              Sfpart;
                typedef typename Ptraits::Degree                        Degree;
                typedef typename Ptraits::Differentiate                 Deriv;
                CGAL::Uncertain<CGAL::Sign> unknown=
                                        Uncertain<CGAL::Sign>::indeterminate();
                CGAL::Uncertain<CGAL::Sign> s=eval_interv(p,
                                                          x.get_left(),
                                                          x.get_right());
                if(!s.is_same(unknown))
                        return s;
                // We are not sure about the sign. We calculate the gcd in
                // order to know if both polynomials have common roots.
                Polynomial_1 sfpp=Sfpart()(p);
                Polynomial_1 gcd=Gcd()(sfpp,Sfpart()(x.get_pol()));
                if(Degree()(gcd)==0)
                        return refine_and_return(p,x);

                // At this point, gcd is not 1; we proceed as follows:
                // -we refine x until having p monotonic in x's interval (to be
                // sure that p has at most one root on that interval),
                // -if the gcd has a root on this interval, both roots are
                // equal (we return 0), otherwise, we refine until having a
                // result.

                // How to assure that p is monotonic in an interval: when its
                // derivative is never zero in that interval.
                Polynomial_1 dsfpp=Deriv()(sfpp);
                CGAL::Gmpfr xl(x.get_left());
                CGAL::Gmpfr xr(x.get_right());
                while(eval_interv(dsfpp,xl,xr).is_same(unknown)){
                        Refiner()(x.get_zpol(),
                                  xl,
                                  xr,
                                  2*CGAL::max(xl.get_precision(),
                                              xr.get_precision()));
                }
                x.set_left(xl);
                x.set_right(xr);

                // How to know that the gcd has a root: evaluate endpoints.
                CGAL::Sign sleft,sright;
                ZSignat sign_at_gcd(PolConverter()(gcd));
                if((sleft=sign_at_gcd(x.get_left()))==CGAL::ZERO||
                   (sright=sign_at_gcd(x.get_right()))==CGAL::ZERO||
                   (sleft!=sright))
                        return CGAL::ZERO;
                return refine_and_return(p,x);
        }
}; // struct Sign_at_z_1

template <class Polynomial_,
          class ZPolynomial_,
          class PolConverter_,
          class Bound_,
          class Algebraic_,
          class Refiner_,
          class Signat_,
          class Ptraits_,
          class ZPtraits_>
class Is_zero_at_z_1:
public CGAL::cpp98::binary_function<Polynomial_,Algebraic_,bool>{
        // This implementation will work with any polynomial type whose
        // coefficient type is explicit interoperable with Gmpfi.
        // TODO: Make this function generic.
        public:
        typedef Polynomial_                                     Polynomial_1;
        typedef ZPolynomial_                                    ZPolynomial_1;
        typedef PolConverter_                                   PolConverter;
        typedef Bound_                                          Bound;
        typedef Algebraic_                                      Algebraic;
        typedef Refiner_                                        Refiner;
        typedef Signat_                                         ZSignat;
        typedef Ptraits_                                        Ptraits;
        typedef ZPtraits_                                       ZPtraits;

        private:
        CGAL::Uncertain<CGAL::Sign> eval_interv(const Polynomial_1 &p,
                                                const Bound &l,
                                                const Bound &r)const{
                typedef typename Ptraits::Substitute                    Subst;
                std::vector<CGAL::Gmpfi> substitutions;
                substitutions.push_back(CGAL::Gmpfi(l,r));
                CGAL::Gmpfi eval=Subst()(p,
                                         substitutions.begin(),
                                         substitutions.end());
                return eval.sign();
        }

        public:
        bool operator()(const Polynomial_1 &p,Algebraic x)const{
                typedef typename Ptraits::Gcd_up_to_constant_factor     Gcd;
                typedef typename Ptraits::Make_square_free              Sfpart;
                typedef typename Ptraits::Degree                        Degree;
                typedef typename Ptraits::Differentiate                 Deriv;
                CGAL::Uncertain<CGAL::Sign> unknown=
                                        Uncertain<CGAL::Sign>::indeterminate();
                CGAL::Uncertain<CGAL::Sign> s=eval_interv(p,
                                                          x.get_left(),
                                                          x.get_right());
                if(!s.is_same(unknown))
                        return (s==CGAL::ZERO);
                // We are not sure about the sign. We calculate the gcd in
                // order to know if both polynomials have common roots.
                Polynomial_1 sfpp=Sfpart()(p);
                Polynomial_1 gcd=Gcd()(sfpp,Sfpart()(x.get_pol()));
                if(Degree()(gcd)==0)
                        return false;

                // At this point, gcd is not 1; we proceed as follows:
                // -we refine x until having p monotonic in x's interval (to be
                // sure that p has at most one root on that interval),
                // -if the gcd has a root on this interval, both roots are
                // equal (we return 0), otherwise, we refine until having a
                // result.

                // How to assure that p is monotonic in an interval: when its
                // derivative is never zero in that interval.
                Polynomial_1 dsfpp=Deriv()(sfpp);
                CGAL::Gmpfr xl(x.get_left());
                CGAL::Gmpfr xr(x.get_right());
                while(eval_interv(dsfpp,xl,xr).is_same(unknown)){
                        Refiner()(x.get_zpol(),
                                  xl,
                                  xr,
                                  2*CGAL::max(xl.get_precision(),
                                              xr.get_precision()));
                }
                x.set_left(xl);
                x.set_right(xr);

                // How to know that the gcd has a root: evaluate endpoints.
                CGAL::Sign sleft,sright;
                ZSignat sign_at_gcd(PolConverter()(gcd));
                return((sleft=sign_at_gcd(x.get_left()))==CGAL::ZERO||
                       (sright=sign_at_gcd(x.get_right()))==CGAL::ZERO||
                       (sleft!=sright));
        }
}; // class Is_zero_at_z_1

// TODO: it says in the manual that this should return a size_type, but test
// programs assume that this is equal to int
template <class Polynomial_,
          class ZPolynomial_,
          class PolConverter_,
          class Isolator_>
struct Number_of_solutions_z_1:
public CGAL::cpp98::unary_function<Polynomial_,int>{
        typedef Polynomial_                                     Polynomial_1;
        typedef ZPolynomial_                                    ZPolynomial_1;
        typedef PolConverter_                                   PolConverter;
        typedef Isolator_                                       Isolator;
        size_t operator()(const Polynomial_1 &p)const{
                // TODO: make sure that p is square free (precondition)
                Isolator isol(PolConverter()(p));
                return isol.number_of_real_roots();
        }
}; // struct Number_of_solutions_z_1

// This functor not only compares two algebraic numbers. In case they are
// different, it refines them until they do not overlap.
template <class Algebraic_,
          class Bound_,
          class Comparator_>
struct Compare_z_1:
public CGAL::cpp98::binary_function<Algebraic_,Algebraic_,CGAL::Comparison_result>{
        typedef Algebraic_                                      Algebraic;
        typedef Bound_                                          Bound;
        typedef Comparator_                                     Comparator;

        CGAL::Comparison_result operator()(Algebraic a,Algebraic b)const{
                Bound al=a.get_left();
                Bound ar=a.get_right();
                Bound bl=b.get_left();
                Bound br=b.get_right();
                CGAL::Comparison_result c=Comparator()(a.get_zpol(),al,ar,
                                                       b.get_zpol(),bl,br);
                a.set_left(al);
                a.set_right(ar);
                b.set_left(bl);
                b.set_right(br);
                return c;
        }

        CGAL::Comparison_result operator()(Algebraic a,const Bound &b)const{
                Bound al=a.get_left();
                Bound ar=a.get_right();
                Algebraic balg(b);
                CGAL::Comparison_result c=Comparator()(a.get_zpol(),al,ar,
                                                       balg.get_zpol(),b,b);
                a.set_left(al);
                a.set_right(ar);
                return c;
        }

        template <class T>
        CGAL::Comparison_result operator()(Algebraic a,const T &b)const{
                Bound al=a.get_left();
                Bound ar=a.get_right();
                Algebraic balg(b);
                CGAL::Comparison_result c=Comparator()(a.get_zpol(),
                                                       al,
                                                       ar,
                                                       balg.get_zpol(),
                                                       balg.get_left(),
                                                       balg.get_right());
                a.set_left(al);
                a.set_right(ar);
                return c;
        }

}; // class Compare_z_1

template <class Algebraic_,
          class Bound_,
          class Comparator_>
struct Bound_between_z_1:
public CGAL::cpp98::binary_function<Algebraic_,Algebraic_,Bound_>{
        typedef Algebraic_                                      Algebraic;
        typedef Bound_                                          Bound;
        typedef Comparator_                                     Comparator;

        Bound operator()(Algebraic a,Algebraic b)const{
                typedef Compare_z_1<Algebraic,Bound,Comparator> Compare;
                typename Bound::Precision_type prec;
                switch(Compare()(a,b)){
                        case CGAL::LARGER:
                                CGAL_assertion(b.get_right()<a.get_left());
                                prec=CGAL::max(b.get_right().get_precision(),
                                               a.get_left().get_precision());
                                        return Bound::add(b.get_right(),
                                                          a.get_left(),
                                                          1+prec)/2;
                                break;
                        case CGAL::SMALLER:
                                CGAL_assertion(a.get_right()<b.get_left());
                                prec=CGAL::max(a.get_right().get_precision(),
                                               b.get_left().get_precision());
                                        return Bound::add(a.get_right(),
                                                          b.get_left(),
                                                          1+prec)/2;
                                break;
                        default:
                                CGAL_error_msg(
                                        "bound between two equal numbers");
                }
        }
}; // struct Bound_between_z_1

template <class Polynomial_,
          class ZPolynomial_,
          class PolConverter_,
          class Bound_,
          class Algebraic_,
          class Isolator_,
          class Comparator_,
          class Signat_,
          class Ptraits_,
          class ZPtraits_>
struct Isolate_z_1:
public CGAL::cpp98::binary_function<Algebraic_,Polynomial_,std::pair<Bound_,Bound_> >{
        typedef Polynomial_                                     Polynomial_1;
        typedef ZPolynomial_                                    ZPolynomial_1;
        typedef PolConverter_                                   PolConverter;
        typedef Bound_                                          Bound;
        typedef Algebraic_                                      Algebraic;
        typedef Isolator_                                       Isolator;
        typedef Comparator_                                     Comparator;
        typedef Signat_                                         Signat;
        typedef Ptraits_                                        Ptraits;
        typedef ZPtraits_                                       ZPtraits;

        std::pair<Bound,Bound>
        operator()(const Algebraic &a,const Polynomial_1 &p)const{
                std::vector<Algebraic> roots;
                std::back_insert_iterator<std::vector<Algebraic> > rit(roots);
                typedef Solve_z_1<Polynomial_1,
                                  ZPolynomial_1,
                                  PolConverter,
                                  Bound,
                                  Algebraic,
                                  Isolator,
                                  Signat,
                                  Ptraits,
                                  ZPtraits>                     Solve;
                typedef Compare_z_1<Algebraic,Bound,Comparator> Compare;
                Solve()(p,false,rit);
                for(typename std::vector<Algebraic>::size_type i=0;
                    i<roots.size();
                    ++i){
                        // we use the comparison functor, that makes both
                        // intervals disjoint iff the algebraic numbers they
                        // represent are not equal
                        Compare()(a,roots[i]);
                }
                return std::make_pair(a.left(),a.right());
        }
}; // Isolate_z_1

template <class Polynomial_,
          class Bound_,
          class Algebraic_,
          class Refiner_>
struct Approximate_absolute_z_1:
public CGAL::cpp98::binary_function<Algebraic_,int,std::pair<Bound_,Bound_> >{
        typedef Polynomial_                                     Polynomial_1;
        typedef Bound_                                          Bound;
        typedef Algebraic_                                      Algebraic;
        typedef Refiner_                                        Refiner;

        // This implementation assumes that Bound is Gmpfr.
        // TODO: make generic.
        std::pair<Bound,Bound> operator()(const Algebraic &x,const int &a)const{
                Bound xl(x.get_left()),xr(x.get_right());
                // refsteps=log2(xl-xr)
                mpfr_t temp;
                mpfr_init(temp);
                mpfr_sub(temp,xr.fr(),xl.fr(),GMP_RNDU);
                mpfr_log2(temp,temp,GMP_RNDU);
                long refsteps=mpfr_get_si(temp,GMP_RNDU);
                mpfr_clear(temp);
                Refiner()(x.get_zpol(),xl,xr,CGAL::abs(refsteps+a));
                x.set_left(xl);
                x.set_right(xr);
                CGAL_assertion(a>0?
                               (xr-xl)*CGAL::ipower(Bound(2),a)<=Bound(1):
                               (xr-xl)<=CGAL::ipower(Bound(2),-a));
                return std::make_pair(xl,xr);
        }
}; // struct Approximate_absolute_z_1

template <class Polynomial_,
          class Bound_,
          class Algebraic_,
          class Refiner_>
struct Approximate_relative_z_1:
public CGAL::cpp98::binary_function<Algebraic_,int,std::pair<Bound_,Bound_> >{
        typedef Polynomial_                                     Polynomial_1;
        typedef Bound_                                          Bound;
        typedef Algebraic_                                      Algebraic;
        typedef Refiner_                                        Refiner;

        std::pair<Bound,Bound> operator()(const Algebraic &x,const int &a)const{
                if(CGAL::is_zero(x))
                        return std::make_pair(Bound(0),Bound(0));
                Bound error=CGAL::ipower(Bound(2),CGAL::abs(a));
                Bound xl(x.get_left()),xr(x.get_right());
                Bound max_b=(CGAL::max)(CGAL::abs(xr),CGAL::abs(xl));
                while(a>0?(xr-xl)*error>max_b:(xr-xl)>error*max_b){
                        Refiner()(x.get_zpol(),
                                  xl,
                                  xr,
                                  std::max<unsigned>(
                                        CGAL::abs(a),
                                        CGAL::max(xl.get_precision(),
                                                  xr.get_precision())));
                        max_b=(CGAL::max)(CGAL::abs(xr),CGAL::abs(xl));
                }
                x.set_left(xl);
                x.set_right(xr);
                CGAL_assertion(
                                a>0?
                                (xr-xl)*CGAL::ipower(Bound(2),a)<=max_b:
                                (xr-xl)<=CGAL::ipower(Bound(2),-a)*max_b);
                return std::make_pair(xl,xr);
        }
}; // struct Approximate_relative_z_1

} // namespace RS_AK1
} // namespace CGAL

#endif // CGAL_RS_FUNCTORS_Z_1_H

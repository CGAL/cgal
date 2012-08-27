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

#ifndef CGAL_RS_FUNCTORS_H
#define CGAL_RS_FUNCTORS_H

#include <CGAL/basic.h>
#include <mpfi.h>
#include <CGAL/RS/polynomial_1.h>
#include <CGAL/RS/algebraic_1.h>
#include <CGAL/RS/polynomial_1_utils.h>
#include <CGAL/RS/solve_1.h>
#include <CGAL/RS/sign_1.h>
#include <CGAL/RS/sign_1_rs.h>
#include <CGAL/RS/refine_1_rs.h>
#include <CGAL/RS/compare_1.h>
#include <CGAL/RS/polynomial_converter.h>
#include <CGAL/Gmpfr.h>

#ifdef IEEE_DBL_MANT_DIG
#  define CGAL_RS_FUNCTORS_DBL_PREC IEEE_DBL_MANT_DIG
#else
#  define CGAL_RS_FUNCTORS_DBL_PREC 53
#endif

namespace RSFunctors{

typedef CGAL::RS_polynomial_1           Polynomial;
typedef CGAL::Algebraic_1               Algebraic;
typedef CGAL::Gmpfr                     Bound;
typedef int                             Multiplicity;

template <class _P>
struct Compute_polynomial_1:
public std::unary_function<Algebraic,_P>{
        typedef _P                      P;
        typedef CGAL::from_rs_poly<P>   back;
        P operator()(const Algebraic &a)const{
                return back()(a.pol());
        }
};      // Compute_polynomial_1

template <class _P,class _Gcd_policy>
struct Is_square_free_1:
public std::unary_function<_P,bool>{
        typedef _P                      P;
        typedef CGAL::to_rs_poly<P>     convert;
        typedef _Gcd_policy             Gcd;
        bool operator()(const P &p)const{
                Polynomial rsp=convert()(p);
                return(!(Gcd()(rsp,rsp.derive()).get_degree_static()));
        }
};      // Is_square_free_1

template <class _P,class _Gcd_policy>
struct Make_square_free_1:
public std::unary_function<_P,_P>{
        typedef _P                      P;
        typedef CGAL::to_rs_poly<P>     convert;
        typedef CGAL::from_rs_poly<P>   back;
        typedef _Gcd_policy             Gcd;
        P operator()(const P &p)const{
                return back()(CGAL::sfpart_1<Gcd>()(convert()(p)));
        }
};      // Make_square_free_1

template <class _P,class _Gcd_policy>
struct Square_free_factorize_1{
        typedef _P                      P;
        typedef CGAL::to_rs_poly<P>     convert;
        typedef CGAL::from_rs_poly<P>   back;
        typedef _Gcd_policy             Gcd;
        template <class OutputIterator>
        OutputIterator operator()(const P &p,OutputIterator oi)const{
                Polynomial rsp=convert()(p);
                CGAL::sqfrvec factorization(CGAL::sqfr_1<Gcd>()(rsp));
                for(CGAL::sqfrvec::iterator i=factorization.begin();
                    i!=factorization.end();
                    ++i){
                        *oi++=std::make_pair(back()((*i).first),(*i).second);
                    }
                return oi;
        }
};      // Square_free_factorize_1

template <class _P,class _Gcd_policy>
struct Is_coprime_1:
public std::binary_function<_P,_P,bool>{
        typedef _P                      P;
        typedef CGAL::to_rs_poly<P>     convert;
        typedef _Gcd_policy             Gcd;
        bool operator()(const P &p1,const P &p2)const{
                CGAL::RS_polynomial_1 rsp1=convert()(p1);
                CGAL::RS_polynomial_1 rsp2=convert()(p2);
                return(!Gcd()(rsp1,rsp2).get_degree_static());
        }
};      // Is_coprime_1

template <class _P,class _Gcd_policy>
struct Make_coprime_1{
        typedef _P                      P;
        typedef CGAL::to_rs_poly<P>     convert;
        typedef CGAL::from_rs_poly<P>   back;
        typedef _Gcd_policy             Gcd;
        bool operator()(const P &p1,const P &p2,P &g,P &q1,P &q2)const{
                typedef _Gcd_policy     Gcd;
                CGAL::RS_polynomial_1 rsp1=convert()(p1);
                CGAL::RS_polynomial_1 rsp2=convert()(p2);
                CGAL::RS_polynomial_1 rsg=convert()(g);
                rsg=Gcd()(rsp1,rsp2);
                g=back()(rsg);
                // even when g==1, we calculate q1 and q2
                q1=back()(*CGAL::Ediv_1()(rsp1,rsg));
                q2=back()(*CGAL::Ediv_1()(rsp2,rsg));
                return rsg.get_degree_static()?false:true;
        }
};      // Make_coprime_1

template <class _Gcd_policy>
struct Solve_RS_1{
    typedef _Gcd_policy         Gcd;
    typedef CGAL::sfpart_1<Gcd> Sfpart;
    typedef CGAL::sqfr_1<Gcd>   Sqfr;

    template <class OutputIterator>
    OutputIterator operator()(const Polynomial &p,OutputIterator res)const{
        int nr,*m;
        mpfi_ptr *x;
        CGAL::sqfrvec sfv=Sqfr()(p);
        x=(mpfi_ptr*)malloc(Sfpart()(p).get_degree()*sizeof(mpfi_ptr));
        m=(int*)calloc(Sfpart()(p).get_degree(),sizeof(int));
        nr=solve_1(x,Sfpart()(p));
        CGAL_assertion_msg(nr>=0,"error in resolution");
        for(CGAL::sqfrvec::size_type i=0;i<sfv.size();++i){
            int k=sfv[i].first.get_degree();
            for(int j=0;k&&j<nr;++j){
                if(!m[j]){
                    CGAL::Sign sg_l=
                        CGAL::RSSign::signat(sfv[i].first,&(x[j]->left));
                    CGAL::Sign sg_r=
                        CGAL::RSSign::signat(sfv[i].first,&(x[j]->right));
                    if((sg_l!=sg_r)||((sg_l==CGAL::ZERO)&&(sg_r==CGAL::ZERO))){
                        m[j]=sfv[i].second;
                        --k;
                    }
                }
            }
        }
        for(int i=0;i<nr;++i){
            *res++=std::make_pair(*new Algebraic(x[i],p,i,m[i]
                                                 //,i?x[i-1]:NULL,
                                                 //i==nr-1?NULL:x[i+1]
                                                 ),
                                  m[i]);
        }
        free(m);
        free(x);
        return res;
    }

    template <class OutputIterator>
    OutputIterator operator()(const Polynomial &p,
                              bool known_to_be_square_free,
                              OutputIterator res)const{
        int nr,m;
        mpfi_ptr *x;
        if(known_to_be_square_free){
            p.set_sf();
            x=(mpfi_ptr*)malloc(p.get_degree()*sizeof(mpfi_ptr));
            nr=solve_1(x,p);
            CGAL_assertion_msg(nr>=0,"error in resolution");
            m=1;    // we know in this case that multiplicity is 1
        }else{
            x=(mpfi_ptr*)malloc(Sfpart()(p).get_degree()*sizeof(mpfi_ptr));
            nr=solve_1(x,Sfpart()(p));
            CGAL_assertion_msg(nr>=0,"error in resolution");
            m=0;    // we won't calculate multiplicities
        }
        for(int i=0;i<nr;++i)
            *res++=*new Algebraic(x[i],p,i,m
                            //,i?x[i-1]:NULL,
                            //i==nr-1?NULL:x[i+1]
                            );
        free(x);
        return res;
    }
};  // Solve_RS_1

template <class _P,class _Gcd_policy>
struct Solve_1{
    typedef _P                  P;
    typedef CGAL::to_rs_poly<P> convert;
    typedef _Gcd_policy         Gcd;
    typedef Solve_RS_1<Gcd>     Solve_RS;

    template <class OutputIterator>
    OutputIterator operator()(const P &p,OutputIterator res)const{
        return Solve_RS()(convert()(p),res);
    }

    template <class OutputIterator>
    OutputIterator operator()(const P &p,
                              bool known_to_be_square_free,
                              OutputIterator res)const{
        return Solve_RS()(convert()(p),known_to_be_square_free,res);
    }

  template <class OutputIterator>
  OutputIterator operator()(const P &p,
                            const Bound& lower,
                            const Bound& upper,
                            OutputIterator res)const{
    typedef std::vector<std::pair<Algebraic,Multiplicity> > RMV;
    RMV roots;
    this->operator()(p,std::back_inserter(roots));
    for(typename RMV::iterator it = roots.begin(); it != roots.end();it++){
      if(lower <= it->first && it->first <= upper)
        *res++=*it;
    }
    return res;
  }

  template <class OutputIterator>
  OutputIterator operator()(const P &p,
                            bool known_to_be_square_free,
                            const Bound& lower,
                            const Bound& upper,
                            OutputIterator res)const{
    typedef std::vector< Algebraic > RV;
    RV roots;
    this->operator()(p,known_to_be_square_free,std::back_inserter(roots));
    for(typename RV::iterator it = roots.begin(); it != roots.end();it++){
      if(lower <= *it && *it <= upper)
        *res++=*it;
    }
    return res;
  }
};  // Solve_1

template <class _P,class _Coeff_type,class _Gcd>
struct Construct_alg_1{
    typedef _P                          Poly;
    typedef _Coeff_type                 Coeff;
    typedef _Gcd                        Gcd;
    typedef CGAL::to_rs_poly<Poly>      convert;
    typedef Solve_1<Poly,Gcd>           Solve;

    Algebraic operator()(int a) const {
        return Algebraic(a);
    }

    Algebraic operator()(const Bound a) const {
        return Algebraic(a);
    }

    Algebraic operator()(const Coeff a) const {
        return Algebraic(a);
    }

    Algebraic operator()(const Poly &p,int i) const {
      CGAL_precondition(CGAL::is_square_free(p));
      std::vector<Algebraic> roots;
      std::back_insert_iterator<std::vector<Algebraic> > rootsit(roots);
      Solve()(p,true,rootsit);
      return roots[i];
    }

    Algebraic operator()(const Poly &p,Bound l,Bound u) const {
        mpfi_t i;
        mpfi_init(i);
        mpfr_set(&i->left,l.fr(),GMP_RNDD);
        mpfr_set(&i->right,u.fr(),GMP_RNDU);
        return Algebraic(i,convert()(p),0,0
                         //,NULL,NULL
                        );
    }
};  // Construct_alg_1

template <class _P>
struct Number_of_solutions_1:
public std::unary_function<_P,int>{
    typedef _P                  P;
    typedef CGAL::to_rs_poly<P> convert;

    int operator()(const P &p)const{
        int nr;
        mpfi_ptr *x;
        CGAL::RS_polynomial_1 rspoly=convert()(p);
        x=(mpfi_ptr*)malloc(rspoly.get_degree()*sizeof(mpfi_ptr));
        nr=solve_1(x,rspoly);
        CGAL_assertion_msg(nr>=0,"error in resolution");
        free(x);
        return nr;
    }
};  // Number_of_solutions_1

template <class _P,class _Gcd_policy>
struct Sign_at_1:
public std::binary_function<_P,Algebraic,CGAL::Sign>{
    typedef _P                  P;
    typedef _Gcd_policy         Gcd;
    typedef CGAL::to_rs_poly<P> convert;

    CGAL::Sign operator()(const P &p,const Algebraic &a)const{
        return RS3::sign_1(convert()(p),a);
    }
};  // Sign_at_1

template <class _P,class _Gcd_policy>
struct Is_zero_at_1:
public std::binary_function<_P,Algebraic,bool>{
    typedef _P                  P;
    typedef _Gcd_policy         Gcd;
    typedef Sign_at_1<P,Gcd>    Sign_at;

    bool operator()(const P &p,const Algebraic &a)const{
        return (Sign_at()(p,a)==CGAL::ZERO);
    }
};  // Is_zero_at_1

template <class _Gcd_policy>
struct Compare_1:
    public std::binary_function<Algebraic,Algebraic,CGAL::Comparison_result>{
  typedef _Gcd_policy                   Gcd;
  typedef CGAL::Comparison_result       Comparison_result;
  typedef CGAL::Gmpz                    Gmpz;
  typedef CGAL::Gmpq                    Gmpq;

  Comparison_result operator()(const Algebraic &r1,const Algebraic &r2)const{
    return CGAL::RS_COMPARE::compare_1<Gcd>(r1,r2);
  }

  Comparison_result operator()(const int &r1,  const Algebraic &r2)const{
    return this->operator()(Algebraic(r1),r2);}
  Comparison_result operator()(const Bound &r1,const Algebraic &r2)const{
    return this->operator()(Algebraic(r1),r2);}
  Comparison_result operator()(const Gmpz &r1, const Algebraic &r2)const{
    return this->operator()(Algebraic(r1),r2);}
  Comparison_result operator()(const Gmpq &r1, const Algebraic &r2)const{
    return this->operator()(Algebraic(r1),r2);}
  Comparison_result operator()(const Algebraic &r1,const int   &r2)const{
    return this->operator()(r1,Algebraic(r2));}
  Comparison_result operator()(const Algebraic &r1,const Bound &r2)const{
    return this->operator()(r1,Algebraic(r2));}
  Comparison_result operator()(const Algebraic &r1,const Gmpz  &r2)const{
    return this->operator()(r1,Algebraic(r2));}
  Comparison_result operator()(const Algebraic &r1,const Gmpq  &r2)const{
    return this->operator()(r1,Algebraic(r2));}
};  // Compare_1

template <class _P,class _Gcd>
struct Isolate_1:
public std::binary_function<Algebraic,_P,std::pair<Bound,Bound> >{
    typedef _P                          Poly;
    typedef _Gcd                        Gcd;
    typedef CGAL::to_rs_poly<Poly>      convert;
    typedef Solve_1<Poly,Gcd>           Solve;
    typedef Compare_1<Gcd>              Compare;

    std::pair<Bound,Bound> operator()(const Algebraic &a,const Poly &p)const{
        std::vector<Algebraic> roots;
        std::back_insert_iterator<std::vector<Algebraic> > rootsit(roots);
        Solve()(p,true,rootsit);
        for(std::vector<Algebraic>::size_type i=0;i<roots.size();++i)
                Compare()(a,roots[i]);
        return std::make_pair(Bound(a.left()),Bound(a.right()));
    }
};  // Isolate_1

template <class _Gcd_policy>
struct Bound_between_1:
    public std::binary_function<Algebraic,Algebraic,Bound>{
        typedef _Gcd_policy     Gcd;
        Bound operator()(const Algebraic &x1,const Algebraic &x2)const{
            double l,r,m;
            switch(CGAL::RS_COMPARE::compare_1<Gcd>(x1,x2)){
                case CGAL::LARGER:
                    CGAL_assertion(x2.sup()<x1.inf());
                    l=x2.sup().to_double(std::round_toward_infinity);
                    r=x1.inf().to_double(std::round_toward_neg_infinity);
                    m=(l+r)/2;
                    if(l<m&&m<r){
                        return Bound(m,CGAL_RS_FUNCTORS_DBL_PREC);
                    }
                    return Bound::add(x2.sup(),
                                      x1.inf(),
                                      (x2.sup().get_precision()>
                                                x1.inf().get_precision()?
                                       1+x2.sup().get_precision():
                                       1+x1.inf().get_precision()))/2;
                    break;
                case CGAL::SMALLER:
                    CGAL_assertion(x1.sup()<x2.inf());
                    l=x1.sup().to_double(std::round_toward_infinity);
                    r=x2.inf().to_double(std::round_toward_neg_infinity);
                    m=(l+r)/2;
                    if(l<m&&m<r){
                        return Bound(m,CGAL_RS_FUNCTORS_DBL_PREC);
                    }
                    return Bound::add(x1.sup(),
                                      x2.inf(),
                                      (x1.sup().get_precision()>
                                                x2.inf().get_precision()?
                                       1+x1.sup().get_precision():
                                       1+x2.inf().get_precision()))/2;
                    break;
                default:
                    CGAL_error_msg("bound between two equal numbers");
            }
        }
    };  // Bound_between_1

struct Approximate_absolute_1:
    public std::binary_function<Algebraic,int,std::pair<Bound,Bound> >{

  std::pair<Bound,Bound>
  operator()(const Algebraic& x, int prec) const {
    RS3::refine_1(x,CGAL::abs(prec));

    CGAL_postcondition_code(
     CGAL::Gmpfr::Precision_type
        subprec=1+
            std::max<CGAL::Gmpfr::Precision_type>(x.inf().get_precision(),
                                                  x.sup().get_precision());
        CGAL::Gmpfr width=CGAL::Gmpfr::sub(x.sup(),x.inf(),subprec);
    )
    CGAL_postcondition_code(if(prec>0))
    CGAL_postcondition(width*CGAL::ipower(Bound(2),prec)<=Bound(1));
    CGAL_postcondition_code(else)
    CGAL_postcondition(width<=CGAL::ipower(Bound(2),-prec));

    return std::make_pair(x.inf(),x.sup());
  }
};

struct Approximate_relative_1
  :public std::binary_function<Algebraic,int,std::pair<Bound,Bound> >{

  std::pair<Bound,Bound> operator()(const Algebraic &x, int prec) const {
    if(CGAL::is_zero(x))
      return std::make_pair(Bound(0),Bound(0));

    Bound error=CGAL::ipower(Bound(2),CGAL::abs(prec));
    Bound max_b=(CGAL::max)(CGAL::abs(x.sup()),CGAL::abs(x.inf()));
    while(prec>0?
          (x.sup()-x.inf())*error>max_b:
          (x.sup()-x.inf())>error*max_b){
      RS3::refine_1(x,mpfi_get_prec(x.mpfi())+CGAL::abs(prec));
      max_b=(CGAL::max)(CGAL::abs(x.sup()),CGAL::abs(x.inf()));
    }

    CGAL_postcondition_code(if(prec>0))
    CGAL_postcondition((x.sup()-x.inf())*CGAL::ipower(Bound(2),prec)<=max_b);
    CGAL_postcondition_code(else)
    CGAL_postcondition((x.sup()-x.inf())<=CGAL::ipower(Bound(2),-prec)*max_b);

    return std::make_pair(x.inf(),x.sup());
  }
};

} // namespace RSFunctors

#undef CGAL_RS_FUNCTORS_DBL_PREC

#endif  // CGAL_RS_FUNCTORS_H

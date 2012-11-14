// Copyright (c) 2007-2009 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_POLYNOMIAL_1_UTILS_H
#define CGAL_RS_POLYNOMIAL_1_UTILS_H

#include <CGAL/basic.h>
#include <gmp.h>
#include <CGAL/RS/polynomial_1.h>
#include <CGAL/RS/polynomial_converter.h>
#include <CGAL/RS/solve_1.h>
#ifdef CGAL_RS_USE_UGCD
#include <CGAL/RS/ugcd/ugcd.h>
#endif
#include <rs_exports.h>
#ifdef CGAL_USE_RS3
#include <rs3_fncts.h>
#endif

namespace CGAL{

#ifdef CGAL_USE_RS3
// Fabrice's modular gcd algorithm
struct Rsgcd_1:
public std::binary_function<RS_polynomial_1,RS_polynomial_1,RS_polynomial_1>{
    RS_polynomial_1& operator()(
                                const RS_polynomial_1 &p1,
                                const RS_polynomial_1 &p2)const{
        int dr,d1,d2;
        mpz_t * r_z;
        d1=p1.get_degree();
        d2=p2.get_degree();
        dr=rs3_up_mz_gcd(
                        &r_z,
                        (const mpz_t*)p1.get_coefs(),
                        d1,
                        (const mpz_t*)p2.get_coefs(),
                        d2);
        RS_polynomial_1 *result=new RS_polynomial_1(&r_z,dr);
        return *result;
    }
};
#endif

struct Cgalgcd_1:
public std::binary_function<RS_polynomial_1,RS_polynomial_1,RS_polynomial_1>{
    RS_polynomial_1& operator()(
                                const RS_polynomial_1 &p1,
                                const RS_polynomial_1 &p2)const{
        typedef Polynomial<Gmpz>    P;
        typedef from_rs_poly<P>     convert;
        typedef to_rs_poly<P>       convertback;
        return convertback()(CGAL::gcd(convert()(p1),convert()(p2)));
    }
};

#ifdef CGAL_RS_USE_UGCD
// my modular gcd algorithm
struct Modgcd_1:
public std::binary_function<RS_polynomial_1,RS_polynomial_1,RS_polynomial_1>{
    RS_polynomial_1& operator()(
                                const RS_polynomial_1 &u,
                                const RS_polynomial_1 &v)const{
        RS_polynomial_1 *result=new RS_polynomial_1(
                u.get_degree()<v.get_degree()?
                u.get_degree_static():
                v.get_degree_static());
        int dr=
            RS_MGCD::Ugcd::ugcd(
                result->get_coefs(),
                u.get_coefs(),
                u.get_degree_static(),
                v.get_coefs(),
                v.get_degree_static());
        result->force_degree(dr);
        return *result;
    }
};
#endif // CGAL_RS_USE_UGCD

// Cont()(c,u) stores in c the gcd of the coefficients of u
struct Cont:
public std::binary_function<mpz_ptr,RS_polynomial_1,void>{
    void operator()(mpz_ptr c,const RS_polynomial_1 &u){
        mpz_set(c,u.coef(0));
        for(int i=1;i<u.get_degree_static()+1;++i)
            mpz_gcd(c,c,u.coef(i));
    }
};

// Pp()(u)=u/Cont()(u)
struct Pp:
public std::unary_function<RS_polynomial_1,RS_polynomial_1>{
    RS_polynomial_1& operator()(const RS_polynomial_1 &u){
        mpz_t c;
        mpz_init(c);
        Cont()(c,u);
        RS_polynomial_1 *res=new RS_polynomial_1(u/c);
        mpz_clear(c);
        return *res;
    }
};

// pseudo division; returns <q,r>, where:
// d^l*f=qg+r, d=leadingcoef(f), l=max(deg(f)-deg(g)+1,0)
struct Pdiv_1:
public std::binary_function<
                            RS_polynomial_1,
                            RS_polynomial_1,
                            std::pair<RS_polynomial_1,RS_polynomial_1> >{
    std::pair<RS_polynomial_1,RS_polynomial_1>
    operator()(const RS_polynomial_1 &f,const RS_polynomial_1 &g){
        int degf,degg,lambda,i;
        mpz_t d;
        mpz_ptr lcg;
        degf=f.get_degree();
        degg=g.get_degree();
        RS_polynomial_1 q(degf);
        RS_polynomial_1 r(f);
        lcg=g.leading_coefficient();
        lambda=degf-degg+1;
        if(lambda>0){
            mpz_init(d);
            mpz_pow_ui(d,g.leading_coefficient(),lambda);
            r*=d;
            mpz_clear(d);
        }
        if(!degg){
            for(int i=0;i<=f.get_degree_static();++i)
                mpz_divexact(q.coef(i),r.coef(i),lcg);
            return std::make_pair(q,RS_polynomial_1());
        }
        // don't use get_degree_static
        while((i=r.get_degree()-degg)>=0){
            mpz_divexact(q.coef(i),r.leading_coefficient(),lcg);
            r-=g.times_monomial(q.coef(i),i);
        }
        return std::make_pair(q,r);
    }
};

struct Pdivrem_1:
public std::binary_function<RS_polynomial_1,RS_polynomial_1,RS_polynomial_1>{
    RS_polynomial_1
    operator()(const RS_polynomial_1 &f,const RS_polynomial_1 &g){
        return Pdiv_1()(f,g).second;
    }
};

struct Pdivquo_1:
public std::binary_function<RS_polynomial_1,RS_polynomial_1,RS_polynomial_1>{
    RS_polynomial_1
    operator()(const RS_polynomial_1 &f,const RS_polynomial_1 &g){
        return Pdiv_1()(f,g).first;
    }
};

// subresultant GCD algorithm:
// Knuth, TAOCP vol. 2, 2nd edition, section 4.6.1, algorithm C
struct Gcd_1:
public std::binary_function<RS_polynomial_1,RS_polynomial_1,RS_polynomial_1>{
    RS_polynomial_1& operator()(
                                const RS_polynomial_1 &uu,
                                const RS_polynomial_1 &vv)const{
        RS_polynomial_1 u,r;
        RS_polynomial_1 *v=new RS_polynomial_1(vv.get_degree());
        mpz_t contu,g,d,h,den_h;
        int delta;
        // C1.
        mpz_init(contu);
        Cont()(contu,uu);
        u=uu/contu;
        mpz_init(g);
        Cont()(g,vv); // g will store cont(v)
        *v=vv/g;
        mpz_init(d);
        mpz_gcd(d,contu,g);
        // we don't need cont(u) and cont(v) anymore; so we will use the
        // variables which stored them (we will use contu as a temp and g to
        // store the value g)
        mpz_set_ui(g,1);    // now, g will store the g from the algorithm
        mpz_init_set_ui(h,1);
        mpz_init(den_h);
        // C2.
        while((r=Pdivrem_1()(u,*v)).get_degree_static()){
            delta=u.get_degree()-v->get_degree();
            // C3.
            u=*v;
            if(delta<0){
                // v := r(x) * (gh^(-delta))
                mpz_pow_ui(contu,h,(unsigned)(-delta));
                mpz_mul(contu,g,contu);
                *v=r*contu;
            }else{
                // v := r(x) / (gh^delta)
                mpz_pow_ui(contu,h,(unsigned)delta);
                mpz_mul(contu,g,contu);
                *v=r/contu;
            }
            mpz_set(g,u.leading_coefficient());
            if(delta){  // if delta=0, h remains unchanged
                if(delta==1){   // if delta=1, h:=g
                    mpz_set(h,g);
                }else{  // otherwise, h:=h^(1-delta)*g^delta
                    if(delta>1){
                        mpz_pow_ui(den_h,h,(unsigned)(delta-1));
                        mpz_pow_ui(h,g,(unsigned)delta);
                        mpz_divexact(h,h,den_h);
                    }else{  // delta<0
                        mpz_pow_ui(h,h,(unsigned)(1-delta));
                        mpz_pow_ui(den_h,g,(unsigned)(-delta));
                        mpz_divexact(h,h,den_h);
                    }
                }
            }
        }
        mpz_clear(contu);
        mpz_clear(g);
        mpz_clear(h);
        mpz_clear(den_h);
        if(mpz_sgn(r.coef(0))){
            v->force_degree(0);
            // v(x)=1, so d*Pp()(v(x))=d
            v->set_coef(0,d);
        }else{
            *v=Pp()(*v);
            *v*=d;
        }
        mpz_clear(d);
        return *v;
    }
};

// exact division: the result is undefined when the division is not exact
struct Ediv_1:
public std::binary_function<RS_polynomial_1,RS_polynomial_1,RS_polynomial_1*>{
    RS_polynomial_1*
    operator()(const RS_polynomial_1 &f,const RS_polynomial_1 &g){
        /*
        int degf,degg,i;
        mpz_ptr lcg;
        mpz_t r;
        degf=f.get_degree();
        degg=g.get_degree();
        RS_polynomial_1 *q=new RS_polynomial_1(degf-degg);
        lcg=g.leading_coefficient();
        if(!degg){
            for(int i=0;i<=degf;++i)
                mpz_divexact(q->coef(i),f.coef(i),lcg);
            return q;
        }
        mpz_init(r);
        mpz_set(r,f.leading_coefficient());
        std::cout<<"f="<<f<<std::endl;
        std::cout<<"g="<<g<<std::endl;
        for(i=degf-degg;i>0;--i){
            //std::cout<<"\ni="<<i;
            mpz_divexact(q->coef(i),r/ *f.coef(i+degg)* /,lcg);
            //--------------------------------------------------
            // mpz_mul(r,q->coef(i),lcg);
            // mpz_sub(r,f.coef(i+degf-degg-1),r);
            //--------------------------------------------------
            mpz_mul(r,q->coef(i),g.coef(degg-1));
            mpz_sub(r,f.coef(i+degf-degg-1),r);
            std::cout<<"i="<<i<<", q[i]="<<q->coef(i)<<", r="<<Gmpz(r)<<std::endl;
        }
        mpz_divexact(q->coef(0),r/ *f.coef(degg)* /,lcg);
        mpz_clear(r);
        //std::cout<<"\nediv("<<f<<","<<g<<") = "<<(*q)<<std::endl;
        return q;
        */
    // naive implementation
        RS_polynomial_1 *ret=new RS_polynomial_1(Pdivquo_1()(f,g));
        return ret;
    }
};

struct Constantpoly_1:
public std::unary_function<long,RS_polynomial_1>{
    RS_polynomial_1&
    operator()(long i){
        RS_polynomial_1 *cpoly=new RS_polynomial_1(0);
        cpoly->set_coef_si(0,i);
        return *cpoly;
    }
};

// squarefree factorization, Yun's algorithm (1976)
// if P=sum(P_i^i), this function returns a vector with pairs (P_i,i)
// it gives for free the square-free part of P, which is C_1
template<class _Gcd_policy>
struct do_sqfr_1:
public std::unary_function<RS_polynomial_1,sqfrvec*>{
    typedef _Gcd_policy         Gcd;
    sqfrvec* operator()(const RS_polynomial_1 &P)const{
        sqfrvec *res=new sqfrvec();
        if(!P.get_degree()){
            res->push_back(
                std::make_pair(Constantpoly_1()(mpz_sgn(P.coef(0))?1:0),1));
            return res;
        }
        RS_polynomial_1 dP(P.derive());
        RS_polynomial_1 G(Pp()(Gcd()(P,dP)));
        if(!G.get_degree()){
            res->push_back(std::make_pair(P,1));
            return res;
        }
        if(!P.has_sfpart()){
            polyptr pp(Ediv_1()(P,G));
            P.set_sfpart(pp);
        }
        RS_polynomial_1 C_i(P.sfpart()); // C_1 is P/G, i.e. the sf part
        RS_polynomial_1 D_i(*Ediv_1()(dP,G)-(C_i.derive()));//D_1=dP/G-dC_1
        RS_polynomial_1 P_i;
        for(int i=1;;++i){
            if(D_i.get_degree_static()||mpz_sgn(D_i.coef(0)))
                P_i=Pp()(Gcd()(C_i,D_i));
            else
                P_i=C_i;
            if(P_i.get_degree_static())
                res->push_back(std::make_pair(P_i,i));
            C_i=*Ediv_1()(C_i,P_i); // C_{i+1}=C_i/P_i
            if(!C_i.get_degree_static())
                return res;
            D_i=*Ediv_1()(D_i,P_i)-(C_i.derive()); // D_{i+1}=D_i/P_i-dC_{i+1}
        }
    }
};

template<class _Gcd_policy>
struct sqfr_1:
public std::unary_function<RS_polynomial_1,sqfrvec>{
    typedef _Gcd_policy         Gcd;
    sqfrvec operator()(const RS_polynomial_1 &P){
        if(!P.has_sqfr()){
            sqfrptr sp(do_sqfr_1<Gcd>()(P));
            P.set_sqfr(sp);
        }
        return P.sqfr();
    }
};

// the square-free part of P is P/gcd(P,dP), but it can be also calculated
// as the product of its sf factors
template<class _Gcd_policy>
struct do_sfpart_1:
public std::unary_function<RS_polynomial_1,RS_polynomial_1*>{
    typedef _Gcd_policy         Gcd;
    RS_polynomial_1* operator()(const RS_polynomial_1 &P)const{
        if(P.get_degree()){
            // TODO: not optimal
            return &Pp()(*Ediv_1()(P,Pp()(Gcd()(P,P.derive()))));
        }else
            return (new RS_polynomial_1(
                            Constantpoly_1()(mpz_sgn(P.coef(0))?1:0)));
    }
};

#ifdef CGAL_USE_RS3
template<>
struct do_sfpart_1<Rsgcd_1>:
public std::unary_function<RS_polynomial_1,RS_polynomial_1*>{
    RS_polynomial_1* operator()(const RS_polynomial_1 &P)const{
        if(P.get_degree()){
            int d_sfp;
            mpz_t* sfp_z;
            d_sfp=rs3_up_mz_sfp(
                               &sfp_z,
                               (const mpz_t*)P.get_coefs(),
                               P.get_degree());
            RS_polynomial_1 *result=new RS_polynomial_1(&sfp_z,d_sfp);
            return result;
        }else
            return (new RS_polynomial_1(
                            Constantpoly_1()(mpz_sgn(P.coef(0))?1:0)));
    }
};
#endif

template<class _Gcd_policy>
struct sfpart_1:
public std::unary_function<RS_polynomial_1,RS_polynomial_1*>{
    typedef _Gcd_policy         Gcd;
    const RS_polynomial_1& operator()(const RS_polynomial_1 &P){
        if(!P.has_sfpart()){
            polyptr pp(do_sfpart_1<Gcd>()(P));
            P.set_sfpart(pp);
        }
        return P.sfpart();
    }
};

} // namespace CGAL

#endif  // CGAL_RS_POLYNOMIAL_1_UTILS_H

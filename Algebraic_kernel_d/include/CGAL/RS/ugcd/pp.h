// Copyright (c) 2007 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS__PP_H
#define CGAL_RS__PP_H

#include <gmp.h>
#include "p.h"
#include "pagealloc.h"
#include <cstdio>

namespace CGAL{
namespace RS_MGCD{

class Prime_polynomial:public Prime,public Page_alloc{
    protected:
        static
        int pp_from_poly(CGALRS_PN *pp,mpz_t *poly,int n){
            int i;
            if(!(pp[n]=mpz_fdiv_ui(poly[n],p_prime())))
                return -1;
            for(i=0;i<n;++i)
                pp[i]=mpz_fdiv_ui(poly[i],p_prime());
            return n;
        };

        static
        void pp_out_str(FILE *stream,CGALRS_PN* pp,int n){
            int i;
            for(i=n;i;--i)
                fprintf(stream,"%d*x^%d+",pp[i],i);
            fprintf(stream,"%d",pp[0]);
            return;
        };

        // Knuth 2; m>=n
        static
        int pp_pdivrem(CGALRS_PN *r,CGALRS_PN *u,int m,CGALRS_PN *v,int n){
            int k,j;
            for(k=0;k<=m;++k)
                r[k]=u[k];
            // division p. 402
            //for(k=m-n;k>=0;--k)
            //  for(j=n+k-1;j>=k;--j)
            //      r[j]=p_sub(r[j],p_mul(qk,v[j-k]));
            // pseudo-division, p. 407
            for(k=m-n;k>=0;--k)
                for(j=n+k-1;j>=0;--j)
                    //r[j]=p_sub(
                    //              p_mul(v[n],r[j]),
                    //              p_mul(r[n+k],(j<k?0:v[j-k])));
                    r[j]=p_submuls(v[n],r[j],r[n+k],(j<k?0:v[j-k]));
            --n;
            while(!r[n]&&n)
                --n;
            return n;
        };

        static
        CGALRS_PN pp_pp(CGALRS_PN *pp,CGALRS_PN *p,int dp){
            int i;
            CGALRS_PN inv,cont=p[dp];
            for(i=0;i<dp;++i)
                cont=p_gcd(cont,p[i]);
            inv=p_inv(cont);
            for(i=0;i<=dp;++i)
                pp[i]=p_mul(p[i],inv);
            return cont;
        };

        // GCL, page 280; da>=db
        static
        int pp_gcd(CGALRS_PN *g,CGALRS_PN *a,int da,CGALRS_PN *b,int db){
            CGALRS_PN *r0,*r1,*r2;
            int i,d0,d1,d2;
            d0=da;
            r0=(CGALRS_PN*)palloc((1+da)*sizeof(CGALRS_PN));
            //for(i=0;i<=da;++i)
            //  r0[i]=a[i];
            pp_pp(r0,a,da);
            d1=db;
            r1=(CGALRS_PN*)palloc((1+da)*sizeof(CGALRS_PN));
            //for(i=0;i<=db;++i)
            //  r1[i]=b[i];
            pp_pp(r1,b,db);
            r2=(CGALRS_PN*)palloc((1+da)*sizeof(CGALRS_PN));
            d2=pp_pdivrem(r2,r0,d0,r1,d1);
            while(d2){
                for(i=0;i<=d1;++i)
                    r0[i]=r1[i];
                d0=d1;
                //for(i=0;i<=d2;++i)
                //  r1[i]=r2[i];
            pp_pp(r1,r2,d2);
            d1=d2;
            d2=pp_pdivrem(r2,r0,d0,r1,d1);
            }
            if(!r2[0]){
                //CGALRS_PN inv=p_inv(r1[d1]);
                //g[d1]=1;
                //for(i=0;i<d1;++i)
                //    g[i]=p_mul(inv,r1[i]);
                for(i=0;i<=d1;++i)
                    g[i]=r1[i];
            }else{
                d1=0;
                g[0]=1;
            }
            CGALRS_PFREE(r0);
            CGALRS_PFREE(r1);
            CGALRS_PFREE(r2);
            return d1;
        };

}; // class Prime_polynomial

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__PP_H

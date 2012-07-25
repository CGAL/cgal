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

#ifndef CGAL_RS__CRT_H
#define CGAL_RS__CRT_H

#include "pp.h"
#include <gmp.h>
#include <vector>
#include <boost/multi_array.hpp>

namespace CGAL{
namespace RS_MGCD{

class Crt:public Prime_polynomial{
    protected:

        // chinese remainder torture (GCL, page 180)
        static
        void cra(mpz_ptr r,CGALRS_PN * m,CGALRS_PN *u,int n){
            int k,i;
            CGALRS_PN product,temp;
            std::vector<CGALRS_PN> v(n);

            v[0]=u[0];
            for(k=1;k<n;++k){
                p_set_prime(m[k]);
                // step 1, gamma_k is p_inv(product)
                product=p_convert(m[0]);
                for(i=1;i<k;++i)
                    //product=p_mul(product,p_convert(m[i]));
                    product=p_mulc(product,m[i]);
                // step 2
                temp=p_convert(v[k-1]);
                for(i=k-2;i>=0;--i)
                    //temp=p_add(p_convert(v[i]),p_mul(temp,p_convert(m[i])));
                    temp=p_mulcaddc(temp,m[i],v[i]);
                //v[k]=p_mul(p_sub(p_convert(u[k]),temp),p_inv(product));
                v[k]=p_convsubdiv(u[k],temp,product);
            }
            // step 3: operations are done in Zm, not in Z
            CGALRS_mpz_set_spn(r,p_pntospn(v[n-1]));
            for(k=n-2;k>=0;--k){
                CGALRS_mpz_mul_pn(r,r,m[k]);
                CGALRS_mpz_add_pn(r,r,v[k]);
            }
            return;
        };

        // polynomial chinese remainder algorithm, it is the same:
        // m are the modules, and m has size size_y, p is the residue
        // vector and also has size size_y, but every one of its elements
        // has size size_x (which will be de degree of the output
        // polynomial);
        // size_y is what is called n in the book
        static
        void pcra(mpz_t *r,
                  CGALRS_PN *m,
                  std::vector<CGALRS_PN* > p,
                  int size_x,
                  int size_y){
            typedef boost::multi_array<CGALRS_PN,2> pn_matrix;
            typedef pn_matrix::index pn_matrix_index;
            pn_matrix v(boost::extents[size_x+1][size_y]);
            pn_matrix_index i,j,k;
            CGALRS_PN product,temp;

            for(j=0;j<=size_x;++j){
                v[j][0]=p[0][j];
                for(k=1;k<size_y;++k){
                    p_set_prime(m[k]);
                    // step 1, gamma_k is p_inv(product)
                    product=p_convert(m[0]);
                    for(i=1;i<k;++i)
                        product=p_mulc(product,m[i]);
                    // step 2
                    temp=p_convert(v[j][k-1]);
                    for(i=k-2;i>=0;--i)
                        temp=p_mulcaddc(temp,m[i],v[j][i]);
                    v[j][k]=p_convsubdiv(p[k][j],temp,product);
                }
            }
            // step 3
            // be careful: operations are done in Zm, not in Z
            for(j=0;j<=size_x;++j){
                CGALRS_mpz_set_spn(r[j],p_pntospn(v[j][size_y-1]));
                for(k=size_y-2;k>=0;--k){
                    CGALRS_mpz_mul_pn(r[j],r[j],m[k]);
                    CGALRS_mpz_add_pn(r[j],r[j],v[j][k]);
                }
            }
        };

}; // class Crt

} // namespace RS_MGCD
} // namespace CGAL

#endif  // CGAL_RS__CRT_H

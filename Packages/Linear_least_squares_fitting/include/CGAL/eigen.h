// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Bruno Levy, Pierre Alliez

#ifndef CGAL_EIGEN_H
#define CGAL_EIGEN_H

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class FT>
void eigen_semi_definite_symmetric(const FT *mat, 
                                   int n, 
                                   FT *eigen_vec, 
                                   FT *eigen_val,
                                   const int MAX_ITER = 100) 
{
  static const FT EPS = 0.00001;
  static const FT PPI = 3.14159265358;
      
  // number of entries in mat
  int nn = (n*(n+1))/2;
      
  // copy matrix
  FT *a = new FT[nn];
  for(int ij=0; ij<nn; ij++) 
    a[ij] = mat[ij];
  // Fortran-porting
  a--;
      
  // init diagonalization matrix as the unit matrix
  FT *v = new FT[n*n];
  ij = 0;
  for(int i=0; i<n; i++)
    for(j=0; j<n; j++) 
      if(i==j)
        v[ij++] = 1.0;
      else
        v[ij++] = 0.0;
  // Fortran-porting
  v--;
      
  // compute weight of the non diagonal terms 
  ij = 1;
  FT a_norm = 0.0;
  for(i=1; i<=n; i++)
    for(j=1; j<=i; j++) 
    {
      if( i!=j ) 
      {
        a_ij = a[ij];
        a_norm += a_ij * a_ij;
      }
      ij++;
    }
      
  if(a_norm != 0.0) 
  {
    FT a_normEPS = a_norm * EPS;
    FT thr = a_norm;
  
    // rotations
    int nb_iter = 0;
    while(thr > a_normEPS && nb_iter < MAX_ITER) 
    {
      nb_iter++;
      FT thr_nn = thr / nn;
          
      for(int l=1; l< n; l++) 
      {
        for(int m=l+1; m<=n; m++) 
        {
          // compute sinx and cosx 
          lq = (l*l-l)/2;
          mq = (m*m-m)/2;
          
          lm = l + mq;
          a_lm = a[lm];
          a_lm_2 = a_lm * a_lm;
          
          if(a_lm_2 < thr_nn)
            continue;
          
          ll   = l+lq;
          mm   = m+mq;
          a_ll = a[ll];
          a_mm = a[mm];
          
          FT delta = a_ll - a_mm;
          
          if(delta == 0.0)
            x = - PPI/4 ; 
          else 
            x = - atan( (a_lm+a_lm) / delta ) / 2.0;

          sinx    = sin(x);
          cosx    = cos(x);
          sinx_2  = sinx * sinx;
          cosx_2  = cosx * cosx;
          sincos  = sinx * cosx;
          
          // rotate L and M columns 
          ilv = n*(l-1);
          imv = n*(m-1);
          
          for( i=1; i<=n;i++ ) 
          {
            if( (i!=l) && (i!=m) ) 
            {
              iq = (i*i-i)/2;
              
              if( i<m )  
                im = i + mq; 
              else
                im = m + iq;
              a_im = a[im];
              
              if( i<l ) 
                il = i + lq; 
              else 
                il = l + iq;
              a_il = a[il];
              
              a[il] = a_il * cosx - a_im * sinx;
              a[im] = a_il * sinx + a_im * cosx;
            }
            
            ilv++;
            imv++;
            
            v_ilv = v[ilv];
            v_imv = v[imv];
            
            v[ilv] = cosx * v_ilv - sinx * v_imv;
            v[imv] = sinx * v_ilv + cosx * v_imv;
          } 
          
          x = a_lm * sincos; 
          x += x;
          
          a[ll] =  a_ll * cosx_2 + a_mm * sinx_2 - x;
          a[mm] =  a_ll * sinx_2 + a_mm * cosx_2 + x;
          a[lm] =  0.0;
          
          thr = fabs(thr - a_lm_2);
        }
      }
    }         
  }
      
  // convert indices and copy eigen values 
  a++;
  for(i=0; i<n; i++) 
  {
    k = i + (i*(i+1))/2;
    eigen_val[i] = a[k];
  }
  delete [] a;
      
  // sort eigen values and vectors 
  index = new int[n];
  for(int i=0; i<n; i++)
    index[i] = i;
      
  for(int i=0; i<(n-1); i++)
  {
    FT x = eigen_val[i];
    int k = i;
        
    for(j=i+1; j<n; j++) 
      if(x < eigen_val[j]) 
      {
        k = j;
        x = eigen_val[j];
      }
        
    eigen_val[k] = eigen_val[i];
    eigen_val[i] = x;
      
    int jj = index[k];
    index[k] = index[i];
    index[i] = jj;
  }


  // save eigen vectors 
  v++; // back to C++
  ij = 0;
  for(int k=0; k<n; k++ ) 
  {
    int ik = index[k]*n;
    for(int i=0; i<n; i++) 
      eigen_vec[ij++] = v[ik++];
  }
  
  delete [] v;
  delete [] index;
}

} // end namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_EIGEN_H

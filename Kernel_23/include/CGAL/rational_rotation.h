// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_RATIONAL_ROTATION_H
#define CGAL_RATIONAL_ROTATION_H

#include <algorithm>
#include <CGAL/number_type_basic.h>

namespace CGAL {

template < class NT >
void
rational_rotation_approximation( const NT &  dirx,     // dir.x()
                                 const NT &  diry,     // dir.y()
                                       NT &  sin_num,  // return
                                       NT &  cos_num,  // return
                                       NT &  denom,    // return
                                 const NT &  eps_num,  // quality_bound
                                 const NT &  eps_den )
{

  const NT& n   = eps_num;
  const NT& d   = eps_den;
  const NT  NT0 = NT(0)  ;
  const NT  NT1 = NT(1)  ;
  CGAL_kernel_precondition( (dirx != NT0) ||  (diry != NT0));
  CGAL_kernel_precondition( n > NT0 );
  CGAL_kernel_precondition( d > NT0 );
  NT & sin = sin_num;
  NT & cos = cos_num;
  NT & den = denom;
  NT   dx = CGAL_NTS abs(dirx);
  NT   dy = CGAL_NTS abs(diry);
  NT   sq_hypotenuse = dx*dx + dy*dy;
  NT   common_part;
  NT   diff_part;
  NT   rhs;
  bool lower_ok;
  bool upper_ok;

  if (dy > dx)
  {
     std::swap (dx,dy);
  }
  // approximate sin = dy / sqrt(sq_hypotenuse)
  // if ( dy / sqrt(sq_hypotenuse) < n/d )
  if (dy * dy * d * d < sq_hypotenuse * n * n)
  {
      cos = NT1;
      sin = NT0;
      den = NT1;
  }
  else
  {
      NT  p;
      NT  q;
      NT  p0 = NT0;
      NT  q0 = NT1;
      NT  p1 = NT1;
      NT  q1 = NT1;

      for(;;)
      {
          p = p0 + p1;
          q = q0 + q1;
          sin = NT(2)*p*q;
          den = p*p + q*q;

      // sanity check for approximation
      //        sin/den < dy/sqrt(hypotenuse) + n/d
      //    &&  sin/den > dy/sqrt(hypotenuse) - n/d
      // ===    sin/den - n/d  <   dy/sqrt(sq_hypotenuse)
      //    &&  sin/den + n/d  >   dy/sqrt(sq_hypotenuse)
      // ===    (sin^2 d^2 + n^2 den^2)sq_hypotenuse - 2... < dy^2 d^2 den^2
      //    &&  (sin^2 d^2 + n^2 den^2)sq_hypotenuse + 2... > dy^2 d^2 den^2

          common_part = (sin*sin*d*d + n*n*den*den)*sq_hypotenuse;
          diff_part   = NT(2)*n*sin*d*den*sq_hypotenuse;
          rhs         = dy*dy*d*d*den*den;

          upper_ok    = (common_part - diff_part < rhs);
          lower_ok    = (common_part + diff_part > rhs);

          if ( lower_ok && upper_ok )
          {
             // if ( (p*p)%2 + (q*q)%2 > NT1)
             // {
             //     sin = p*q;
             //     cos = (q*q - p*p)/2;    // exact division
             //     den = (p*p + q*q)/2;    // exact division
             // }
             // else
             // {
                    cos = q*q - p*p;
             // }

             break;
          }
          else
          {
              // if ( dy/sqrt(sq_hypotenuse) < sin/den )
              if ( dy*dy*den*den < sin*sin*sq_hypotenuse )
              {
                  p1 = p;
                  q1 = q;
              }
              else
              {
                  p0 = p;
                  q0 = q;
              }
          }
      } // for(;;)
  }
  dx = dirx;
  dy = diry;

  if (CGAL_NTS abs(dy) > CGAL_NTS abs(dx) ) { std::swap (sin,cos); }

  if (dx < NT0) { cos = - cos; }

  if (dy < NT0) { sin = - sin; }

  sin_num = sin;
  cos_num = cos;
  denom   = den;
}


template < class NT >
void
rational_rotation_approximation( const double& angle,
                                            NT &  sin_num,  // return
                                            NT &  cos_num,  // return
                                            NT &  denom,    // return
                                      const NT &  eps_num,  // quality_bound
                                      const NT &  eps_den )
{

  const NT& n   = eps_num;
  const NT& d   = eps_den;
  const NT  NT0 = NT(0)  ;
  const NT  NT1 = NT(1)  ;
  CGAL_kernel_precondition( n > NT0 );
  CGAL_kernel_precondition( d > NT0 );
  NT& isin = sin_num;
  NT& icos = cos_num;
  NT& iden = denom;
  double dsin = std::sin(angle);
  double dcos = std::cos(angle);
  double dn = CGAL::to_double(n);
  double dd = CGAL::to_double(d);
  double eps = dn / dd;
  dsin = CGAL_NTS abs( dsin);
  dcos = CGAL_NTS abs( dcos);
  NT   common_part;
  NT   diff_part;
  NT   os;
  bool lower_ok;
  bool upper_ok;
  bool swapped = false;

  if (dsin > dcos)
  {
     swapped = true;
     std::swap (dsin,dcos);
  }
  if ( dsin < eps )
  {
      icos = NT1;
      isin = NT0;
      iden = NT1;
  }
  else
  {
      NT  p;
      NT  q;
      NT  p0 = NT0;
      NT  q0 = NT1;
      NT  p1 = NT1;
      NT  q1 = NT1;

      for(;;)
      {
          p = p0 + p1;
          q = q0 + q1;
          isin = NT(2)*p*q;
          iden = p*p + q*q;

          // XXX sanity check for approximation
          //        sin/den < dsin + n/d
          //    &&  sin/den > dsin - n/d
          //        sin < dsin * den + n/d * den
          //    &&  sin > dsin * den - n/d * den
          os          = CGAL::to_double(isin);
          diff_part   = eps  * CGAL::to_double(iden);
          common_part = dsin * CGAL::to_double(iden);

          upper_ok    = (common_part - diff_part < os);
          lower_ok    = (os < common_part + diff_part);

          if ( lower_ok && upper_ok )
          {
             // if ( (p*p)%2 + (q*q)%2 > NT1)
             // {
             //     isin = p*q;
             //     icos = (q*q - p*p)/2;    // exact division
             //     iden = (p*p + q*q)/2;    // exact division
             // }
             // else
             // {
                    icos = q*q - p*p;
             // }

             break;
          }
          else
          {
              // XXX if ( dsin < sin/den )
              if ( dsin * CGAL::to_double(iden) < CGAL::to_double(isin) )
              {
                  p1 = p;
                  q1 = q;
              }
              else
              {
                  p0 = p;
                  q0 = q;
              }
          }
      } // for(;;)
  }

  if ( swapped ) { std::swap (isin,icos); }

  dsin = std::sin( angle);
  dcos = std::cos( angle);
  if (dcos < 0.0) { icos = - icos; }
  if (dsin < 0.0) { isin = - isin; }

  sin_num = isin;
  cos_num = icos;
  denom   = iden;
}

} //namespace CGAL

#endif // CGAL_RATIONAL_ROTATION_H

/* Copyright 2004
   Stanford University

   This file is part of the DSR PDB Library.

   The DSR PDB Library is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

   The DSR PDB Library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the DSR PDB Library; see the file LICENSE.LGPL.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA. */


#ifndef CGAL_DSRPDB_RMS_H
#define CGAL_DSRPDB_RMS_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/Matrix.h>
#include <CGAL/PDB/range.h>
#include <cmath>
#include <vector>
namespace CGAL { namespace PDB {


//! Compute the cRMS of the collection of Points after transforming the first
template <class RangeA, class RangeB>
double cRMS(const RangeA& ra, const RangeB& rb,
            const Transform &tr=Transform(1,0,0,0,
                                          0,1,0,0,
                                          0,0,1,0)) {
  CGAL_precondition(CGAL::PDB::distance(ra) == CGAL::PDB::distance(rb));

  double ret=0;
  int num=0;
  {
    typename RangeB::iterator bc= rb.begin();
    typename RangeA::iterator ac= ra.begin();
    for (; bc != rb.end(); ++bc, ++ac){
      //std::cout << *ac << " " << *bc << std::endl;
      double cd= squared_distance(tr(*ac), *bc);
      //std::cout << cd << std::endl;
      ret+=cd;
      ++num;
    }
  }
  return std::sqrt(ret/ num);
}



//! Compute the dRMS of the collection of Points without alignment
template <class RangeA, class RangeB>
double dRMS(RangeA ra, RangeB rb) {
  CGAL_assertion(ra.size()==rb.size());
  double ret=0;
  int count=0;
  typename RangeA::iterator ac= ra.begin();
  typename RangeB::iterator bc= rb.begin();
  for (; ac != ra.end(); ++ac, ++bc) {
    typename RangeA::iterator ac2= ra.begin();
    typename RangeB::iterator bc2= rb.begin();
    for (; ac2 != ac; ++ac2, ++bc2) {
      double da= std::sqrt((*ac-*ac2).squared_length());
      double db= std::sqrt((*bc-*bc2).squared_length());
      ret+= (da-db)*(da-db);
      ++count;
    }
  }
  return ret/count;
}

//! Return the distance matrix
template <class RangeA, class RangeB>
Matrix distance_matrix(RangeA ra, RangeB rb) {
  int dista= std::distance(ra.begin(), ra.end());
  int distb= std::distance(rb.begin(), rb.end());
  Matrix ret(dista, distb);
  
  int indi=0;
  for (typename RangeA::iterator ac= ra.begin(); ac != ra.end(); ++ac, ++indi){
    int indj=0;
    for (typename RangeB::iterator ac2= rb.begin(); ac2!= rb.end(); ++ac2, ++indj){
      Vector dif= *ac- *ac2;
      //std::cout << dif << std::endl;
      double d= std::sqrt(dif*dif);
      //std::cout << d << std::endl;
      ret[indi][indj]=d;
    }
  }
  return ret;
}

//! Return the distance matrix
template <class Range>
Matrix distance_matrix(Range r) {
    return distance_matrix(r,r);
}



}}
#endif

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
#include <CGAL/PDB/geometry.h>
#include <CGAL/PDB/Matrix.h>
#include <cmath>
#include <vector>
CGAL_PDB_BEGIN_NAMESPACE
//! Compute the cRMS of the collection of Points without alignment
template <class ItA, class ItB>
double cRMS(ItA ba, ItA ea, ItB bb, ItB eb) {
  Squared_distance sd;
  if (std::distance(ab, ae) != std::distance(bb, be)){
    CGAL_PDB_INTERNAL_NS::error_logger.new_fatal_error("Protein chains used for computing cRMS must have equal lengths.\n");
    return std::numeric_traits<double>::infinity();
  }
  double ret=0;
  int num=0;
  for (It bc= bb, ac= ab; bc != be; ++bc, ++ac){
    
    Point pt= Point(*bc);
    Point tpt= f(pt);
    ret += sd(*ac, tpt);
    ++num;
  }
  return std::sqrt(ret)/ num;
}

//! Compute the dRMS of the collection of Points without alignment
template <class ItA, class ItB>
double dRMS(ItA ba, ItA ea, ItB bb, ItB eb) {
  CGAL_assertion(std::distance(ba, ea) == std::distance(bb, eb));
  double ret=0;
  int count=0;
  for (It ac= ba, bc= bb; ac != ea; ++ac, ++bc) {
    for (It ac2= ba, bc2= bb; ac2 != ac; ++ac2, ++bc2) {
      Vector va= *ac- *ac2;
      Vector vb= *bc- *bc2;
      double da= std::sqrt(va*va);
      double db= std::sqrt(vb*vb);
      ret+= (da-db)*(da-db);
      ++count;
    }
  }
  return ret/count;
}

//! Return the distance matrix
template <class ItA, class ItB>
Matrix distance_matrix(ItA ba, ItA ea, ItB bb, ItB eb) {
  int dist= std::distance(ab, ae);
  Matrix ret(std::distance(ba, ea),std::distance(bb, e));
  
  int indi=0;
  for (ItA ac= ba; ac != ea; ++ac, ++indi){
    int indj=0;
    for (ItB ac2= bb; ac2!= eb; ++ac2, ++indj){
      Vector dif= *ac- *ac2;
      //std::cout << dif << std::endl;
      double d= std::sqrt(dif*dif);
      //std::cout << d << std::endl;
      ret[indi][indj]=d;
    }
  }
  return ret;
}

CGAL_PDB_END_NAMESPACE
#endif

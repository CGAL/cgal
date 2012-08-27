// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_DESCARTES_ROOT_COUNT_H
#define CGAL_POLYNOMIAL_DESCARTES_ROOT_COUNT_H

#include <CGAL/Polynomial/basic.h>

/*! \file
  This file has the root counter used by the Descartes (and other) solvers.
*/

namespace CGAL { namespace POLYNOMIAL { namespace internal {;

class Descartes_root_count
{
protected:
public:
  typedef enum DRC {INVALID=-12, UNKNOWN =-3,ZERO=0, ONE=1, EVEN=2, ODD=3, SOME=5}
  Count;
  Descartes_root_count(unsigned int i) {
    // catch overflow
    CGAL_Polynomial_exactness_precondition(i<10000);
    /*if (i==0) ct_= ZERO;
      else if (i==1) ct_= ONE;
      else if (i%2 == 0) ct_= EVEN;
      else ct_=ODD;*/
    ct_=i;
  }
  Descartes_root_count(Count ct): ct_(ct) {

  }
  Descartes_root_count(): ct_(INVALID){}
  int count() const
  {
    CGAL_precondition(ct_ != INVALID);
    return ct_;
  }
  bool is_unknown() const
  {
    return ct_== UNKNOWN;
  }
  bool is_odd() const
  {
    //CGAL_precondition(ct_ != INVALID);
    if (ct_== INVALID) return false;
    //CGAL_precondition(ct_ != SOME);
    return ct_%2==1;                      // || ct_== SINGLE_ODD;
  }
  template <class O>
  bool operator==(const O &o) const
  {
    bool why_are_you_using_this_function;
    CGAL_precondition(ct_ != INVALID);
    CGAL_precondition(o.ct_ != INVALID);
    return ct_== o.ct_;
  }
  template <class O>
  bool operator!=(const O &o) const
  {
    bool why_are_you_using_this_function;
    CGAL_precondition(ct_ != INVALID);
    CGAL_precondition(o.ct_ != INVALID);
    return ct_ != o.ct_;
  }
  bool is_single() const
  {
    //if (ct_== INVALID) return false;// CGAL_precondition(ct_ != INVALID);
    return ct_== 1;                       // || ct_ == SINGLE_EVEN || ct_ == SINGLE_EVEN;
  }
  bool is_zero() const
  {
    //CGAL_precondition(ct_ != INVALID);
    return ct_== 0;                       // || ct_ == SINGLE_EVEN || ct_ == SINGLE_EVEN;
  }
  static  Descartes_root_count zero() {
    return Descartes_root_count(ZERO);
  }
  static  Descartes_root_count one() {
    return Descartes_root_count(ONE);
  }
  /*static  Descartes_root_count single_even() {
    return Descartes_root_count(SINGLE_EVEN);
    }
    static  Descartes_root_count single_odd() {
    return Descartes_root_count(SINGLE_ODD);
    }*/
  static  Descartes_root_count even() {
    return Descartes_root_count(EVEN);
  }
  static  Descartes_root_count odd() {
    return Descartes_root_count(ODD);
  }
  static  Descartes_root_count some() {
    return Descartes_root_count(SOME);
  }
protected:
  int ct_;
};

inline std::ostream &operator<<(std::ostream &out, Descartes_root_count ct)
{
  switch(ct.count()) {
  case Descartes_root_count::ZERO:
    out <<"ZERO";
    break;
  case Descartes_root_count::ONE:
    out <<"ONE";
    break;
    /*case Descartes_root_count::SINGLE_EVEN:
      out <<"ONE(E)";
      break;
      case Descartes_root_count::SINGLE_ODD:
      out <<"ONE(O)";
      break;*/
  case Descartes_root_count::EVEN:
    out <<"EVEN";
    break;
  case Descartes_root_count::ODD:
    out <<"ODD";
    break;
  case Descartes_root_count::SOME:
    out <<"SOME";
    break;
  default:
    out << "??";
  }
  return out;
}


} } } //namespace CGAL::POLYNOMIAL::internal;
#endif

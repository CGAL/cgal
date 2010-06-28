// Copyright (c) 2009,2010 Inria Lorraine (France). All rights reserved.
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
// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@loria.fr>

#include <CGAL/basic.h>

#if defined(CGAL_USE_GMP) && defined(CGAL_USE_MPFI) && defined(CGAL_USE_RS)

#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>
#include "include/CGAL/_test_algebraic_kernel_1.h"


int main(){
  typedef CGAL::Algebraic_kernel_rs_gmpz_d_1              AK;
  typedef AK::Polynomial_1 Polynomial_1;
  typedef AK::Coefficient Coefficient;
  typedef AK::Bound Bound;
  typedef AK::Algebraic_real_1 Algebraic_real_1;
  typedef AK::Multiplicity_type Multiplicity_type;

  AK ak; // an object of Algebraic_kernel_rs_gmpz_d_1
  CGAL::test_algebraic_kernel_1<AK>(ak);

  AK::Solve_1 solve_1 = ak.solve_1_object();
  Polynomial_1 x = CGAL::shift(AK::Polynomial_1(1),1);
  int returnvalue=0;

  // variant using a bool indicating a square free polynomial
  // multiplicities are not computed
  std::vector<Algebraic_real_1> roots;
  solve_1(x*x-2,true, std::back_inserter(roots));
  if(roots.size()!=2){
    returnvalue-=1;
    std::cerr<<"error 1: the number of roots of x^2-2 must be 2"<<
      std::endl;
  }
  if(-1.42>=roots[0] || -1.41<=roots[0] ||
      1.41>=roots[1] || 1.42<=roots[1]){
    returnvalue-=2;
    std::cerr<<"error 2: the roots of x^2-2 are wrong"<<std::endl;
  }
  roots.clear();

  // variant for roots in a given range of a square free polynomial
  solve_1((x*x-2)*(x*x-3),true, Bound(0),Bound(10),
      std::back_inserter(roots));
  if(roots.size()!=2){
    returnvalue-=4;
    std::cerr<<"error 3: the number of roots of (x^2-2)*(x^2-3)"<<
      " between 0 and 10 must be 2"<<std::endl;
  }
  if(1.41>=roots[0] || 1.42<=roots[0] ||
      1.73>=roots[1] || 1.74<=roots[1]){
    returnvalue-=8;
    std::cerr<<"error 4: the roots of (x^2-2)*(x^2-3)"<<
      " between 0 and 10 are wrong"<<std::endl;
  }
  roots.clear();

  // variant computing all roots with multiplicities
  std::vector<std::pair<Algebraic_real_1,Multiplicity_type> > mroots;
  solve_1((x*x-2), std::back_inserter(mroots));
  if(mroots.size()!=2){
    returnvalue-=16;
    std::cerr<<"error 5: the number of roots of x^2-2 must be 2"<<
      std::endl;
  }
  if(-1.42>=mroots[0].first || -1.41<=mroots[0].first ||
      1.41>=mroots[1].first || 1.42<=mroots[1].first){
    returnvalue-=32;
    std::cerr<<"error 6: the roots of x^2-2 are wrong"<<std::endl;
  }
  if(mroots[0].second!=1 && mroots[1].second!=1){
    returnvalue-=64;
    std::cerr<<"error 7: the multiplicities of the"<<
      " roots of x^2-2 are wrong"<<std::endl;
  }
  mroots.clear();

  // variant computing roots with multiplicities for a range
  solve_1((x*x-2)*(x*x-3),Bound(0),Bound(10),std::back_inserter(mroots));
  if(mroots.size()!=2){
    returnvalue-=128;
    std::cerr<<"error 8: the number of roots of (x^2-2)*(x^2-3)"<<
      " between 0 and 10 must be 2"<<std::endl;
  }
  if(1.41>=mroots[0].first || 1.42<=mroots[0].first ||
      1.73>=mroots[1].first || 1.74<=mroots[1].first){
    returnvalue-=256;
    std::cerr<<"error 9: the roots of (x^2-2)*(x^2-3) are wrong"<<
      std::endl;
  }
  if(mroots[0].second!=1 && mroots[1].second!=1){
    returnvalue-=512;
    std::cerr<<"error 10: the multiplicities of the roots of"<<
      " (x^2-2)*(x^2-3) are wrong"<<std::endl;
  }

  return returnvalue;
}
#else
int main(){
        return 0;
}
#endif

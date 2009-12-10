// Copyright (c) 2009 Inria Lorraine (France). All rights reserved.
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

#ifdef CGAL_USE_RS
#include <CGAL/Algebraic_kernel_d_1_RS_Gmpz.h>

typedef CGAL::Algebraic_kernel_d_1_RS_Gmpz              AK;
typedef AK::Polynomial_1                                Polynomial_1;
typedef AK::Algebraic_real_1                            Algebraic_real_1;
typedef AK::Coefficient                                 Coefficient;
typedef AK::Bound                                       Bound;

#define CGAL_TEST_ALGEBRAIC_REAL_IO(_f) \
        CGAL::set_ascii_mode(ss); \
        CGAL_TEST_ALGEBRAIC_REAL_IO_MODE(_f) \
        CGAL::set_pretty_mode(ss); \
        CGAL_TEST_ALGEBRAIC_REAL_IO_MODE(_f)

#define CGAL_TEST_ALGEBRAIC_REAL_IO_MODE(_f) \
        alg1=_f; \
        ss<<CGAL::oformat(alg1); \
        ss>>CGAL::iformat(alg2); \
        assert(alg1==alg2);

int main(){
  AK ak; // an object of Algebraic_kernel_d_1_RS_Gmpz
  AK::Construct_algebraic_real_1 construct_algreal_1 =
        ak.construct_algebraic_real_1_object();

  Algebraic_real_1 alg1,alg2;
  std::stringstream ss;

  // test construction from int, Coefficient and Bound
  CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algreal_1(int(2)))
  CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algreal_1(Coefficient(2)))
  CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algreal_1(Bound(2)))

  // construction by index
  Polynomial_1 x = CGAL::shift(AK::Polynomial_1(1),1); // the monom x
  CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algreal_1(x*x-2,1))

  // construction by isolating interval
  CGAL_TEST_ALGEBRAIC_REAL_IO(construct_algreal_1(x*x-2,Bound(0),Bound(2)))

  return 0;
}
#else
int main(){
        std::cerr<<"RS was not configured"<<std::endl;
        return 0;
}
#endif

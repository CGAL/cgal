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

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpfr.h>

template<class _NT>
int test_arithmetic(){
        typedef CGAL::Gmpfr     Gmpfr;
        typedef _NT             NT;
        Gmpfr a(7.0);
        NT b(2);
        Gmpfr::set_default_rndmode(std::round_to_nearest);
        if(
              CGAL::Gmpfr::add(a,b) == 9
           && CGAL::Gmpfr::add(a,b,std::round_to_nearest) == 9
           && CGAL::Gmpfr::add(a,b,2) == 8
           && CGAL::Gmpfr::add(a,b,2,std::round_toward_infinity) == 12
           && CGAL::Gmpfr::sub(a,b) == 5
           && CGAL::Gmpfr::sub(a,b,std::round_to_nearest) == 5
           && CGAL::Gmpfr::sub(a,b,2) == 4
           && CGAL::Gmpfr::sub(a,b,2,std::round_toward_infinity) == 6
           && CGAL::Gmpfr::mul(a,b) == 14
           && CGAL::Gmpfr::mul(a,b,std::round_to_nearest) == 14
           && CGAL::Gmpfr::mul(a,b,2) == 16
           && CGAL::Gmpfr::mul(a,b,2,std::round_toward_neg_infinity) == 12
           && CGAL::Gmpfr::div(a,b) == 3.5
           && CGAL::Gmpfr::div(a,b,std::round_to_nearest) == 3.5
           && CGAL::Gmpfr::div(a,b,2) == 4
           && CGAL::Gmpfr::div(a,b,2,std::round_toward_neg_infinity) == 3
           && Gmpfr(9.0).is_square(a) && a==3
        )
                return 0;
        else
                exit(-1);
}

int main(){
        test_arithmetic<CGAL::Gmpfr>();
        test_arithmetic<long>();
        test_arithmetic<unsigned long>();
        test_arithmetic<int>();
        test_arithmetic<CGAL::Gmpz>();
        test_arithmetic<CGAL::Gmpq>();
        return 0;
}

#else
int main(){
        return 0;
}
#endif

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab

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

#ifdef CGAL_USE_MPFI
#include <CGAL/Gmpfi.h>

template<class _NT>
int test_arithmetic(){
        typedef CGAL::Gmpfi     Gmpfi;
        typedef _NT             NT;
        Gmpfi a(7.0);
        NT b(2);
        if(
              CGAL::Gmpfi::add(a,b)            == 9
           && CGAL::Gmpfi::add(b,a,2).inf()    == 8
           && CGAL::Gmpfi::sub(a,b)            == 5
           && CGAL::Gmpfi::sub(b,a,2).sup()    == -4
           && CGAL::Gmpfi::mul(a,b)            == 14
           && CGAL::Gmpfi::mul(b,a,2).sup()    == 16
           && CGAL::Gmpfi::div(a,b)            == 3.5
           && CGAL::Gmpfi::div(b,a,2).inf()    == .25
           && Gmpfi(9.0).is_square(a) && a==3
        )
                return 0;
        else{
                std::cerr<<"error: arithmetic test"<<std::endl;
                exit(-1);
        }
}

int main(){
        test_arithmetic<CGAL::Gmpfi>();
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

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

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpfr.h>
#include <ctime>

template<class _NT>
int test_nt(){
        typedef CGAL::Gmpfr     Gmpfr;
        typedef _NT             NT;
        Gmpfr a(5.0);
        NT b(2);
        if((-a)==-5&&(a+b)==7&&(a-b)==3&&(a*b)==10&&(a/b)==2.5)
                return 0;
        else
                exit(-1);
}

int main(){
        test_nt<CGAL::Gmpfr>();
        test_nt<CGAL::Gmpz>();
        test_nt<CGAL::Gmpq>();
        test_nt<int>();
        test_nt<long>();
        test_nt<unsigned long>();
        return 0;
}

#else
int main(){
        return 0;
}
#endif

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab

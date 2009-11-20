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
#include <CGAL/_test_io.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpzf.h>

int main(){
        typedef CGAL::Gmpfr     Gmpfr;

        Gmpfr plus_infinity;
        Gmpfr minus_infinity;
        CGAL::Gmpzf f;

        mpfr_set_inf(plus_infinity.fr(),1);
        mpfr_set_inf(minus_infinity.fr(),-1);

        CGAL::test_io<Gmpfr,int>(2145338339);
        CGAL::test_io<Gmpfr,double>(.2147483647);
        CGAL::test_io<Gmpfr,Gmpfr>(plus_infinity);
        CGAL::test_io<Gmpfr,Gmpfr>(minus_infinity);
        CGAL::test_io<Gmpfr,CGAL::Gmpzf>(f);

        return 0;
}

#else
int main(){
        return 0;
}
#endif

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab

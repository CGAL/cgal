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
#include <CGAL/_test_io.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>
#include <CGAL/Gmpzf.h>

int main(){
        typedef CGAL::Gmpfr     Gmpfr;
        typedef CGAL::Gmpfi     Gmpfi;

        Gmpfr plus_infinity;
        Gmpfr minus_infinity;
        CGAL::Gmpzf f;

        mpfr_set_inf(plus_infinity.fr(),1);
        mpfr_set_inf(minus_infinity.fr(),-1);

        CGAL::test_interval_io<Gmpfi,int>(2145338339);
        CGAL::test_interval_io<Gmpfi,double>(.2147483647);
        CGAL::test_interval_io<Gmpfi,Gmpfr>(Gmpfr(1,100)/Gmpfr(3,100));
        CGAL::test_interval_io<Gmpfi,std::pair<int,Gmpfr> >
                (std::make_pair(2145338339,plus_infinity));
        CGAL::test_interval_io<Gmpfi,std::pair<Gmpfr,double> >
                (std::make_pair(minus_infinity,.2147483647));
        CGAL::test_interval_io<Gmpfi,std::pair<Gmpfr,Gmpfr> >
                (std::make_pair(minus_infinity,plus_infinity));

        return 0;
}

#else
int main(){
        return 0;
}
#endif

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab

// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Athanasios Kakargias <grad0460@di.uoa.gr>
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (CGAL - Effective Computational Geometry for Curves and Surfaces)

// file : include/CGAL/Root_of/leda_real.h

#ifndef CGAL_ROOT_OF_LEDA_REAL_H
#define CGAL_ROOT_OF_LEDA_REAL_H

#include <CGAL/leda_real.h>
#include <CGAL/Root_of/Root_of_traits.h>

namespace CGAL {

    inline
    leda_real
    make_root_of_2(const leda_real &a, const leda_real &b,
                   const leda_real &c, bool d)
    {
      return CGALi::make_root_of_2_sqrt(a,b,c,d);
    }

    template <>
    struct Root_of_traits< leda_real >
    {
      typedef leda_real RootOf_1;
      typedef leda_real RootOf_2;
      typedef leda_real RootOf_3;
      typedef leda_real RootOf_4;
    };

}

#endif // CGAL_ROOT_OF_LEDA_REAL_H

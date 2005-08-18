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

// file : include/CGAL/Root_of/CGAL_Interval_nt.h

#ifndef CGAL_ROOT_OF_CGAL_INTERVAL_NT_H
#define CGAL_ROOT_OF_CGAL_INTERVAL_NT_H

#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Root_of_2.h>

namespace CGAL {

    // TODO: The CGAL::Interval_nt< P > throws an exception when dealing with 
    // unsafe comparisons that is not catched in the predicates.
  
    template <bool P >
    inline
    CGAL::Interval_nt< P >
    make_root_of_2(const CGAL::Interval_nt<P> &a, const CGAL::Interval_nt<P> &b,
                   const CGAL::Interval_nt<P> &c, bool d)
    {
      return CGALi::make_root_of_2_sqrt(a,b,c,d);
    }

    template <bool P >
    struct Root_of_traits< CGAL::Interval_nt< P > >
    {
      typedef CGAL::Interval_nt< P > RootOf_1;
      typedef CGAL::Interval_nt< P > RootOf_2;
      typedef CGAL::Interval_nt< P > RootOf_3;
      typedef CGAL::Interval_nt< P > RootOf_4;
    };

}

#endif // CGAL_ROOT_OF_CGAL_INTERVAL_NT_H

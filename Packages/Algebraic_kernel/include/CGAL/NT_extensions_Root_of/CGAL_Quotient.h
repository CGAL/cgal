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

// file : include/CGAL/Root_of/CGAL_Quotient.h

#ifndef CGAL_ROOT_OF_CGAL_QUOTIENT_H
#define CGAL_ROOT_OF_CGAL_QUOTIENT_H

#include <CGAL/Quotient.h>
//#include <CGAL/Root_of/Root_of_traits.h>

namespace CGAL {

    template < class NT >
    inline
    typename Root_of_traits< NT >::RootOf_2
    make_root_of_2(const Quotient< NT > &a, const Quotient< NT > &b,
                   const Quotient< NT > &c, bool d)
    {
      return CGALi::make_root_of_2_rational< NT, Quotient< NT > >(a,b,c,d);
    }

    // CGAL::Quotient<NT> should be the same as Root_of_traits<NT>::RootOf_1
    // i.e the default implementation.
    // 
    template < class NT >
    struct Root_of_traits< Quotient< NT > >
      : public Root_of_traits< NT > {};
}

#endif // CGAL_ROOT_OF_CGAL_QUOTIENT_H

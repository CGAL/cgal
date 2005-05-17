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

// file : include/CGAL/Root_of/CGAL_Gmpq.h

#ifndef CGAL_ROOT_OF_CGAL_GMPQ_H
#define CGAL_ROOT_OF_CGAL_GMPQ_H

#include <CGAL/Gmpq.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Gmpz.h>

namespace CGAL {

    inline
    Root_of_traits< CGAL::Gmpz >::RootOf_2
    make_root_of_2(const CGAL::Gmpq &a, const CGAL::Gmpq &b,
                   const CGAL::Gmpq &c, bool d)
    {
      return CGALi::make_root_of_2_rational< CGAL::Gmpz, CGAL::Gmpq >(a,b,c,d);
    }

    template< >
    struct Root_of_traits< CGAL::Gmpq >
      : public Root_of_traits< CGAL::Gmpz > {};
      // CGAL::Gmpq is the same as Root_of_traits< CGAL::Gmpz >::RootOf_1

}

#endif // CGAL_ROOT_OF_CGAL_GMPQ_H

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

// file : include/CGAL/Root_of/CGAL_Gmpz.h

#ifndef CGAL_ROOT_OF_CGAL_GMPZ_H
#define CGAL_ROOT_OF_CGAL_GMPZ_H

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
//#include <CGAL/Root_of/Root_of_traits.h>
#include <CGAL/Root_of_2.h>

namespace CGAL {

    template <>
    struct Root_of_traits< CGAL::Gmpz >
    {
      typedef CGAL::Gmpq RootOf_1;
      typedef Root_of_2< CGAL::Gmpz > RootOf_2;
      typedef Root_of_3< CGAL::Gmpz > RootOf_3;
      typedef Root_of_4< CGAL::Gmpz > RootOf_4;
    };

}

#endif // CGAL_ROOT_OF_CGAL_GMPZ_H

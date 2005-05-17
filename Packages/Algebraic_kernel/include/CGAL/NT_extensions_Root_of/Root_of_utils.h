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

// file : include/CGAL/Root_of/Root_of_utils.h

#ifndef CGAL_ROOT_OF_ROOT_OF_UTILS_H
#define CGAL_ROOT_OF_ROOT_OF_UTILS_H

namespace CGAL{
  namespace CGALi{
   
  // Mini helper to prevent clashes when RT == int (see CGAL/Quotient.h).
//     template < typename T >
//     struct Int_if_not_int { typedef int type; };

//     template <>
//     struct Int_if_not_int<int> { struct type{}; };
    
    #define CGAL_CK_int(T) typename CGALi::Int_if_not_int<T>::type
  }
}

#endif // CGAL_ROOT_OF_ROOT_OF_UTILS_H


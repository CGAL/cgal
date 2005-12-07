// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Athanasios Kakargias
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

// file : include/CGAL/Root_of/double.h

#ifndef CGAL_ROOT_OF_DOUBLE_H
#define CGAL_ROOT_OF_DOUBLE_H


#include <CGAL/Root_of_2.h>


namespace CGAL {

    inline
    double
    make_root_of_2(const double &a, const double &b,
                   const double &c, bool d)
    {
      return CGALi::make_root_of_2_sqrt(a,b,c,d);
    }

    template <>
    struct Root_of_traits< double >
    {
      typedef double RootOf_1;
      typedef double RootOf_2;
      typedef double RootOf_3;
      typedef double RootOf_4;
    };

}

#endif // CGAL_ROOT_OF_DOUBLE_H

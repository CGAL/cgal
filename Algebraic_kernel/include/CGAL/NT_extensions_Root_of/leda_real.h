// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Athanasios Kakargias

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

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

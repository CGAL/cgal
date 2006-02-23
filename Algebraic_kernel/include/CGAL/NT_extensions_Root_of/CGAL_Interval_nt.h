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

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

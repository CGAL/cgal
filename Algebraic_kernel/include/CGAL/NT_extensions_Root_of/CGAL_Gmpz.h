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

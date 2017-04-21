// Copyright (c) 1997   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Sylvain Pion

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/license/Triangulation_2.h>

#ifdef CGAL_NO_DEPRECATED_CODE
#error The class Regular_triangulation_euclidean_traits_2<K> is deprecated;
       the kernel K can be used directly as traits since the weighted point and
       the function objects for weighted points are part of the concept Kernel.
#endif

namespace CGAL {

template < class R, class W = typename R::RT>
class Regular_triangulation_euclidean_traits_2
  : public R
{
public:
  Regular_triangulation_euclidean_traits_2() {}
  Regular_triangulation_euclidean_traits_2(const R& k) : R(k) {}
};

} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

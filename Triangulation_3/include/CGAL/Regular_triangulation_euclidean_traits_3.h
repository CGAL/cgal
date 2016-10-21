// Copyright (c) 1999,2004   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion
//                 Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H

#include <CGAL/basic.h>

namespace CGAL {



template < class K, class Weight = typename K::RT >
class Regular_triangulation_euclidean_traits_3
  : public K
{
  K k;
public:
  Regular_triangulation_euclidean_traits_3(K k) : k(k) {}

  typedef K                                          Kernel;
  typedef typename K::FT                             FT;
  typedef typename K::Point_3                        Bare_point;
  typedef typename K::Weighted_point_3               Weighted_point;
  typedef Weighted_point                             Weighted_point_3;
  typedef Weighted_point                             Point_3;
};


} //namespace CGAL



#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H

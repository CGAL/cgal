// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Stefan Schirra

// This file's name must begin with a lower-case letter for backward 
// compatability.  Unfortunately, you can't have a file that differs only 
// in capitalization on the Windows platforms.

#ifndef CGAL_CONVEX_HULL_TRAITS_2_H
#define CGAL_CONVEX_HULL_TRAITS_2_H

#include <CGAL/license/Convex_hull_2.h>

#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/distance_predicates_2.h>

namespace CGAL {

template <class K_>
class Convex_hull_traits_2 : public K_
{
public:
  Convex_hull_traits_2() { }
  Convex_hull_traits_2(const K_& k) : K_(k) { }
};

template <class K_>
class convex_hull_traits_2 : public Convex_hull_traits_2<K_>
{
  convex_hull_traits_2() { }
  convex_hull_traits_2(const K_& k) : K_(k) { }
};

} //namespace CGAL

#endif // CGAL_CONVEX_HULL_TRAITS_2_H

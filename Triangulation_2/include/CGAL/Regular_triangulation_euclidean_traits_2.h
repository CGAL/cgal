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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Sylvain Pion

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/license/Triangulation_2.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Regular_triangulation_euclidean_traits_2.h>"
#define CGAL_DEPRECATED_MESSAGE_DETAILS \
  "The kernel K can be used directly as traits since weighted points and "\
  "the associated function objects are now part of the concept Kernel."
#include <CGAL/internal/deprecation_warning.h>

namespace CGAL {

template < class K, class W = typename K::RT>
class Regular_triangulation_euclidean_traits_2
  : public K
{
public:
  Regular_triangulation_euclidean_traits_2() {}
  Regular_triangulation_euclidean_traits_2(const K& k) : K(k) {}
};

} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

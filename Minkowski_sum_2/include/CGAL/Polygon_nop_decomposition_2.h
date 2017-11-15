// Copyright (c) 2013  Tel-Aviv University (Israel).
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
// Author(s) : Efi Fogel   <efifogel@gmail.com>

#ifndef CGAL_POLYGON_NOP_DECOMPOSITION_2_H
#define CGAL_POLYGON_NOP_DECOMPOSITION_2_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Polygon_2.h>

namespace CGAL {

/*! \class
 * Nop decomposition strategy.
 * Used for polygons that are already convex.
 */
template <typename Kernel_,
          typename Container_ = std::vector<typename Kernel_::Point_2> >
class Polygon_nop_decomposition_2 {
public:
  typedef Kernel_                                        Kernel;
  typedef Container_                                     Container;
  typedef CGAL::Polygon_2<Kernel, Container>             Polygon_2;

  template <typename OutputIterator_>
  OutputIterator_ operator()(const Polygon_2& pgn, OutputIterator_ oi) const
  { *oi++ = pgn; return oi; }
};

} //namespace CGAL

#endif

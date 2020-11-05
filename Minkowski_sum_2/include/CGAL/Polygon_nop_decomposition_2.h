// Copyright (c) 2013  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_POLYGON_CONVEX_DECOMPOSITION_H
#define CGAL_POLYGON_CONVEX_DECOMPOSITION_H

#include <CGAL/Minkowski_sum_2/Decomposition_strategy_adapter.h>
#include <vector>

namespace CGAL {

/*!
 * \class
 * The O(n^4) optimal strategy for decomposing a polygon into convex
 * sub-polygons.
 */
template <class Kernel_, 
          class Container_ = std::vector<typename Kernel_::Point_2> >
class Optimal_convex_decomposition_2 :
  public Polygon_decomposition_strategy_adapter<Kernel_, Container_,
						Tag_optimal_convex_parition>
{
public:

  typedef Kernel_                                  Kernel;
  typedef CGAL::Polygon_2<Kernel, Container_>      Polygon_2;
  typedef typename Kernel::Point_2                 Point_2;

};

/*!
 * \class
 * Hertel and Mehlhorn's O(n) approximation strategy for decomposing a
 * polygon into convex sub-polygons.
 */
template <class Kernel_, 
          class Container_ = std::vector<typename Kernel_::Point_2> >
class Hertel_Mehlhorn_convex_decomposition_2 :
  public Polygon_decomposition_strategy_adapter<Kernel_, Container_,
						Tag_approx_convex_parition>
{
public:

  typedef Kernel_                                  Kernel;
  typedef CGAL::Polygon_2<Kernel, Container_>      Polygon_2;
  typedef typename Kernel::Point_2                 Point_2;

};

/*!
 * \class
 * Greene's O(n log(n)) approximation strategy for decomposing a polygon into
 * convex sub-polygons.
 */
template <class Kernel_, 
          class Container_ = std::vector<typename Kernel_::Point_2> >
class Greene_convex_decomposition_2 :
  public Polygon_decomposition_strategy_adapter<Kernel_, Container_,
						Tag_Greene_convex_parition>
{
public:

  typedef Kernel_                                  Kernel;
  typedef CGAL::Polygon_2<Kernel, Container_>      Polygon_2;
  typedef typename Kernel::Point_2                 Point_2;

};

} //namespace CGAL

#endif

// Copyright (c) 2006 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein           <wein_r@yahoo.com>

#ifndef CGAL_DECOMPOSITION_STRATEGY_ADAPTER_H
#define CGAL_DECOMPOSITION_STRATEGY_ADAPTER_H

#include <CGAL/partition_2.h>

namespace CGAL {

struct Tag_optimal_convex_parition
{
  bool    dummy;
};

struct Tag_approx_convex_parition
{
  bool    dummy;
};

struct Tag_Greene_convex_parition
{
  bool    dummy;
};

/*!
 * \class
 * An adapter of the global planar polygonal partitioning functions
 * to a decomposition strategy-class.
 */
template <class Kernel_, class Container_, class StrategyTag_>
class Decomposition_strategy_adapter
{
public:

  typedef Kernel_                                  Kernel;
  typedef Polygon_2<Kernel, Container_>            Polygon_2;
  typedef StrategyTag_                             Strategy_tag;

  /*! Default constructor. */
  Decomposition_strategy_adapter ()
  {}

  /*!
   * Decompose a simple polygon to convex sub-polygons.
   * \param pgn The input polygon.
   * \param oi An output iterator of convex polygons.
   * \return A past-the-end iterator for the sub-polygons.
   */
  template <class OutputIterator>
  OutputIterator decompose (const Polygon_2& pgn,
                            OutputIterator oi) const
  {
    // Make a local copy of the polygon, and if it is not counterclockwise
    // oriented, reverse the order of its vertices.
    Polygon_2        my_pgn = pgn;

    if (my_pgn.orientation() == CLOCKWISE)
      my_pgn.reverse_orientation();

    // Perform the decomposition.
    return (_decompose (pgn, Strategy_tag(), oi));
  }

private:

  /*!
   * Decompose the given counter clockwise-oriented polygon using the optimal
   * convex-partition method.
   */
  template <class OutputIterator>
  OutputIterator _decompose (const Polygon_2& pgn,
                             Tag_optimal_convex_parition ,
                             OutputIterator oi) const
  {
    return (optimal_convex_partition_2 (pgn.vertices_begin(),
                                        pgn.vertices_end(),
                                        oi,
                                        Kernel()));
  }

  /*!
   * Decompose the given counter clockwise-oriented polygon using the
   * approximated convex-partition method.
   */
  template <class OutputIterator>
  OutputIterator _decompose (const Polygon_2& pgn,
                             Tag_approx_convex_parition ,
                             OutputIterator oi) const
  {
    return (approx_convex_partition_2 (pgn.vertices_begin(),
                                       pgn.vertices_end(),
                                       oi,
                                       Kernel()));
  }

  /*!
   * Decompose the given counter clockwise-oriented polygon using Greene's
   * approximated convex-partition method.
   */
  template <class OutputIterator>
  OutputIterator _decompose (const Polygon_2& pgn,
                             Tag_Greene_convex_parition ,
                             OutputIterator oi) const
  {
    return (greene_approx_convex_partition_2 (pgn.vertices_begin(),
                                              pgn.vertices_end(),
                                              oi,
                                              Kernel()));
  }
};

} //namespace CGAL

#endif

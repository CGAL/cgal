// Copyright (c) 2014 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ilker O. Yaz

#ifndef CGAL_SURFACE_MODELING_SPOKES_AND_RIMS_ITERATOR_H
#define CGAL_SURFACE_MODELING_SPOKES_AND_RIMS_ITERATOR_H

#include <CGAL/license/Surface_mesh_deformation.h>

/// @cond CGAL_DOCUMENT_INTERNAL

namespace CGAL {
namespace internal {
/**
 * Currently this class is not used by surface modeling package, just leave it for possible future need.
 * Provide simple functionality for iterating over spoke and rim edges
 *   - use get_descriptor() to obtain active edge
 *   - get_iterator() always holds spoke edges */
 /// \code
 /// // how to use Spokes_and_rims_iterator
 /// boost::tie(e_begin, e_end) = out_edges(vertex, halfedge_graph);
 /// Spokes_and_rims_iterator<HalfedgeGraph> rims_it(e_begin, halfedge_graph);
 ///
 /// for ( ; rims_it.get_iterator() != e_end; ++rims_it )
 /// {
 ///   halfedge_descriptor active_hedge = rims_it.get_descriptor();
 ///   // use active_edge as you like
 /// }
 /// \endcode
template<class HalfedgeGraph>
class Spokes_and_rims_iterator
{
public:
  typedef typename boost::graph_traits<HalfedgeGraph>::out_edge_iterator  out_edge_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor  halfedge_descriptor;

  Spokes_and_rims_iterator(out_edge_iterator edge_iterator, HalfedgeGraph& halfedge_graph)
    : is_current_rim(false), iterator(edge_iterator), descriptor(halfedge(*edge_iterator)), halfedge_graph(halfedge_graph)
  { }

  /// descriptor will be assigned to next valid edge, note that iterator might not change
  Spokes_and_rims_iterator<HalfedgeGraph>&
  operator++()
  {
    // loop through one spoke then one rim edge
    if(!is_current_rim && !is_border(descriptor, halfedge_graph)) // it is rim edge's turn
    {
      is_current_rim = true;
      descriptor = next(descriptor, halfedge_graph);
    }
    else // if current edge is rim OR there is no rim edge (current spoke edge is boudary)
    {    // then iterate to next spoke edge
      is_current_rim = false;
      descriptor = halfedge(*(++iterator));
    }
    return *this;
  }

  out_edge_iterator get_iterator()   { return iterator; }
  edge_descriptor   get_descriptor() { return descriptor; }

private:
  bool is_current_rim;        ///< current descriptor is rim or spoke
  out_edge_iterator iterator; ///< holds spoke edges (i.e. descriptor is not always = *iterator)
  halfedge_descriptor descriptor; ///< current active halfedge descriptor for looping
  HalfedgeGraph& halfedge_graph;
};

}//namespace internal
/// @endcond
}//namespace CGAL
#endif //CGAL_SURFACE_MODELING_SPOKES_AND_RIMS_ITERATOR_H

// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/Indirect_edge_compare.h
// package       : $CGAL_Package: Partition_2 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Planar Polygon Partitioning
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Comparison of edges represented by circulators to endpoints. 
// ============================================================================

#ifndef CGAL_INDIRECT_EDGE_COMPARE_H
#define CGAL_INDIRECT_EDGE_COMPARE_H

namespace CGAL {

//
// given circulators to endpoints of two edges, sorts the edges that come
// next (in the direction of circulator) from right to left. This ordering 
// makes finding the edge directly left of a given edge (needed for the 
// y-monotone decomposition algorithm) easy. 
//
template <class Traits>
class Indirect_edge_compare 
{
   public:
     typedef typename Traits::Compare_y_2        Compare_y_2;
     typedef typename Traits::Compare_x_2        Compare_x_2;
     typedef typename Traits::Construct_line_2   Construct_line_2;
     typedef typename Traits::Line_2             Line_2;

     Indirect_edge_compare() : 
          _compare_y_2(Traits().compare_y_2_object()),
          _compare_x_2(Traits().compare_x_2_object()),
          _construct_line_2(Traits().construct_line_2_object())
     { }
     
     // determines if the edge (edge_vtx_1, edge_vtx_1++) has a larger
     // x value than vertex.x() at y-value vertex.y()
     template <class ForwardCirculator>
     bool
     larger_x_at_vertex_y(ForwardCirculator edge_vtx_1, 
                          ForwardCirculator vertex) const;

     template <class ForwardCirculator>
     bool 
     operator()(ForwardCirculator p, ForwardCirculator q) const;

   private:
     Compare_y_2   _compare_y_2;
     Compare_x_2   _compare_x_2;
     Construct_line_2 _construct_line_2;
};

}

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Indirect_edge_compare.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_INDIRECT_EDGE_COMPARE_H


// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Implements the concepts SurfaceMeshCellBase_3, with respect to Mesh_3 spec.
//******************************************************************************

#ifndef MESH_SURFACE_CELL_BASE_3_H_
#define MESH_SURFACE_CELL_BASE_3_H_


#include <CGAL/Regular_triangulation_cell_base_3.h>

namespace CGAL {

namespace Mesh_3 {


/**
 * @class Mesh_surface_cell_base_3
 */
template <class GT,
          class MT,
          class Cb = CGAL::Regular_triangulation_cell_base_3<GT> >
class Mesh_surface_cell_base_3
: public Cb
{
public:
  // Indices
  typedef typename MT::Surface_index  Surface_index;
  typedef typename MT::Index          Index;

  // Triangulation types
  typedef typename Cb::Triangulation_data_structure   Tds;
  typedef typename Tds::Vertex_handle                 Vertex_handle;
  typedef typename Tds::Cell_handle                   Cell_handle;
  typedef typename GT::Point_3                        Point;

  // To get correct cell type in TDS
  template < class TDS3 >
  struct Rebind_TDS
  {
    typedef typename Cb::template Rebind_TDS<TDS3>::Other Cb3;
    typedef Mesh_surface_cell_base_3 <GT, MT, Cb3> Other;
  };


  /// Constructors
  Mesh_surface_cell_base_3()
    : Cb()
    , surface_index_table_()
    , surface_center_table_()
    , bits_(0) { };

  Mesh_surface_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                           Vertex_handle v2, Vertex_handle v3)
    : Cb (v0, v1, v2, v3)
    , surface_index_table_()
    , surface_center_table_()
    , bits_(0) { };

  Mesh_surface_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                           Vertex_handle v2, Vertex_handle v3,
                           Cell_handle n0, Cell_handle n1,
                           Cell_handle n2, Cell_handle n3)
    : Cb (v0, v1, v2, v3, n0, n1, n2, n3)
    , surface_index_table_()
    , surface_center_table_()
    , bits_(0) { };


  /// Destructor
  ~Mesh_surface_cell_base_3() { };

  // Default copy constructor and assignment operator are ok

  /// Set surface index of \c facet to \c index
  void set_surface_index(const int facet, const Surface_index& index)
  {
    CGAL_precondition(facet>=0 && facet<4);
    surface_index_table_[facet] = index;
  }

  /// Returns surface index of facet \c facet
  Surface_index surface_index(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return surface_index_table_[facet];
  }

  /// Marks \c facet as visited
  void set_facet_visited (const int facet)
  {
    CGAL_precondition(facet>=0 && facet <4);
    bits_ |= (1 << facet);
  }

  /// Marks \c facet as not visited
  void reset_visited (const int facet)
  {
    CGAL_precondition(facet>=0 && facet<4);
    bits_ &= (15 & ~(1 << facet));
  }

  /// Returns \c true if \c facet is marked as visited
  bool is_facet_visited (const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return ( (bits_ & (1 << facet)) != 0 );
  }

  /// Sets surface center of \c facet to \c point
  void set_facet_surface_center(const int facet, const Point& point)
  {
    CGAL_precondition(facet>=0 && facet<4);
    surface_center_table_[facet] = point;
  }

  /// Returns surface center of \c facet
  Point get_facet_surface_center(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return surface_center_table_[facet];
  }

  /// Sets surface center index of \c facet to \c index
  void set_facet_surface_center_index(const int facet, const Index& index)
  {
    CGAL_precondition(facet>=0 && facet<4);
    surface_center_index_table_[facet] = index;
  }

  /// Returns surface center of \c facet
  Index get_facet_surface_center_index(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return surface_center_index_table_[facet];
  }

  /// Returns true if facet lies on a surface patch
  bool is_facet_on_surface(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return ( Surface_index() != surface_index_table_[facet]);
  }

private:
  /// Stores surface_index for each facet of the cell
  Surface_index surface_index_table_[4];
  /// Stores surface center of each facet of the cell
  Point surface_center_table_[4];
  /// Stores surface center index of each facet of the cell
  Index surface_center_index_table_[4];
  /// Stores visited facets (4 first bits)
  char bits_;

};  // end class Mesh_surface_cell_base_3




}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // MESH_SURFACE_CELL_BASE_3_H_

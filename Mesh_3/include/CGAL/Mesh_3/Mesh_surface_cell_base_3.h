// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Implements the concepts SurfaceMeshCellBase_3, with respect to Mesh_3 spec.
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_SURFACE_CELL_BASE_3_H
#define CGAL_MESH_3_MESH_SURFACE_CELL_BASE_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Mesh_3/io_signature.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/atomic.h>
#endif

#ifdef _MSC_VER
// Kill warning "C4351: new behavior: elements of array
// 'CGAL::Mesh_3::Mesh_surface_cell_base_3<GT,MT,Cb>::surface_index_table_'
// will be default initialized"
#  pragma warning(disable:4351)
#endif

namespace CGAL {

namespace Mesh_3 {

  
/************************************************
// Class Mesh_surface_cell_base_3_base
// Two versions: sequential / parallel
************************************************/

// Sequential
template <typename Concurrency_tag>
class Mesh_surface_cell_base_3_base
{
public:
  Mesh_surface_cell_base_3_base() 
    : bits_(0) {}
  
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

protected:
  /// Stores visited facets (4 first bits)
  char bits_;
};

#ifdef CGAL_LINKED_WITH_TBB
// Parallel
template<>
class Mesh_surface_cell_base_3_base<Parallel_tag>
{
public:
  Mesh_surface_cell_base_3_base() 
  {
    bits_ = 0;
  }
  
  /// Marks \c facet as visited
  void set_facet_visited (const int facet)
  {
    CGAL_precondition(facet>=0 && facet<4);
    char current_bits = bits_;
    while (bits_.compare_and_swap(current_bits | (1 << facet), current_bits) != current_bits)
    {
      current_bits = bits_;
    }
  }

  /// Marks \c facet as not visited
  void reset_visited (const int facet)
  {
    CGAL_precondition(facet>=0 && facet<4);
    char current_bits = bits_;
    while (bits_.compare_and_swap(current_bits & (15 & ~(1 << facet)), current_bits) != current_bits)
    {
      current_bits = bits_;
    }
  }

  /// Returns \c true if \c facet is marked as visited
  bool is_facet_visited (const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return ( (bits_ & (1 << facet)) != 0 );
  }

protected:
  /// Stores visited facets (4 first bits)
  tbb::atomic<char> bits_;
};
#endif // CGAL_LINKED_WITH_TBB


/**
 * @class Mesh_surface_cell_base_3
 */
template <class GT,
          class MT,
          class Cb>
class Mesh_surface_cell_base_3
: public Mesh_surface_cell_base_3_base<
    typename Cb::Triangulation_data_structure::Concurrency_tag>
, public Cb
{
public:
  // Indices
  typedef typename MT::Surface_patch_index  Surface_patch_index;
  typedef typename MT::Index                Index;

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
  { }

  Mesh_surface_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                           Vertex_handle v2, Vertex_handle v3)
    : Cb (v0, v1, v2, v3)
    , surface_index_table_()
    , surface_center_table_()
  { }

  Mesh_surface_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                           Vertex_handle v2, Vertex_handle v3,
                           Cell_handle n0, Cell_handle n1,
                           Cell_handle n2, Cell_handle n3)
    : Cb (v0, v1, v2, v3, n0, n1, n2, n3)
    , surface_index_table_()
    , surface_center_table_()
  { }


  /// Destructor
  ~Mesh_surface_cell_base_3() { }

  // Default copy constructor and assignment operator are ok

  /// Set surface index of \c facet to \c index
  void set_surface_patch_index(const int facet, const Surface_patch_index& index)
  {
    CGAL_precondition(facet>=0 && facet<4);
    surface_index_table_[facet] = index;
  }

  /// Returns surface index of facet \c facet
  Surface_patch_index surface_patch_index(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return surface_index_table_[facet];
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
    return ( Surface_patch_index() != surface_index_table_[facet]);
  }
  
  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;
  
  void set_surface_index(const int facet, const Surface_index& index)
  { set_surface_patch_index(facet,index); }
  
  /// Returns surface index of facet \c facet
  Surface_index surface_index(const int facet) const
  { return surface_patch_index(facet); }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------

  static
  std::string io_signature()
  {
    return
      Get_io_signature<Cb>()() + "+(" +
      Get_io_signature<Surface_patch_index>()() + ")[4]";
  }
private:
  /// Stores surface_index for each facet of the cell
  Surface_patch_index surface_index_table_[4];
  /// Stores surface center of each facet of the cell
  Point surface_center_table_[4];
  /// Stores surface center index of each facet of the cell
  Index surface_center_index_table_[4];

};  // end class Mesh_surface_cell_base_3

#ifdef _MSC_VER
#  pragma warning(default:4351)
#endif

template < class GT, class MT, class Cb >
inline
std::istream&
operator>>(std::istream &is, Mesh_surface_cell_base_3<GT, MT, Cb> &c)
{
  typename Mesh_surface_cell_base_3<GT, MT, Cb>::Surface_patch_index index;
  is >> static_cast<Cb&>(c);
  for(int i = 0; i < 4; ++i)
  {
    if(is_ascii(is))
      is >> index;
    else
    {
      read(is, index);
    }
    c.set_surface_patch_index(i, index);
  }
  return is;
}

template < class GT, class MT, class Cb >
inline
std::ostream&
operator<<(std::ostream &os,
           const Mesh_surface_cell_base_3<GT, MT, Cb> &c)
{
  os << static_cast<const Cb&>(c);
  for(int i = 0; i < 4; ++i)
  {
    if(is_ascii(os))
      os << ' ' << c.surface_patch_index(i);
    else
      write(os, c.surface_patch_index(i));
  }
  return os;
}


}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_MESH_SURFACE_CELL_BASE_3_H

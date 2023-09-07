// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008,2011 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau, Stephane Tayeb, Andreas Fabri, Jane Tournois

#ifndef CGAL_COMPACT_SIMPLICIAL_MESH_CELL_BASE_3_H
#define CGAL_COMPACT_SIMPLICIAL_MESH_CELL_BASE_3_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/basic.h>
#include <CGAL/Has_timestamp.h>
#include <CGAL/tags.h>
#include <CGAL/TDS_3/internal/Dummy_tds_3.h>

#include <CGAL/SMDS_3/io_signature.h>

namespace CGAL {

// Adds information to Cb about the cell of the input complex containing it
template <typename Subdomain_index_,
          typename Surface_patch_index_,
          typename TDS>
class Compact_simplicial_mesh_cell_3
{
public:
  using Triangulation_data_structure = TDS;

  using Vertex_handle = typename TDS::Vertex_handle;
  using Cell_handle = typename TDS::Cell_handle;
  using Vertex = typename TDS::Vertex;
  using Cell = typename TDS::Cell;
  using TDS_data = typename TDS::Cell_data;

  // Index Type
  using Subdomain_index = Subdomain_index_;
  using Surface_patch_index = Surface_patch_index_;

public:
  // Constructors
  Compact_simplicial_mesh_cell_3()
    : time_stamp_(std::size_t(-1))
  {}

  Compact_simplicial_mesh_cell_3(const Compact_simplicial_mesh_cell_3& rhs)
    : N(rhs.N)
    , V(rhs.V)
    , time_stamp_(rhs.time_stamp_)
    , subdomain_index_(rhs.subdomain_index_)
  {
    for(int i=0; i <4; i++){
      surface_index_table_[i] = rhs.surface_index_table_[i];
    }
  }

  Compact_simplicial_mesh_cell_3(Vertex_handle v0,
                                 Vertex_handle v1,
                                 Vertex_handle v2,
                                 Vertex_handle v3)
    : V(CGAL::make_array(v0, v1, v2, v3))
  {
  }

  Compact_simplicial_mesh_cell_3(Vertex_handle v0,
                                 Vertex_handle v1,
                                 Vertex_handle v2,
                                 Vertex_handle v3,
                                 Cell_handle n0,
                                 Cell_handle n1,
                                 Cell_handle n2,
                                 Cell_handle n3)
    : N(CGAL::make_array(n0, n1, n2, n3))
    , V(CGAL::make_array(v0, v1, v2, v3))
  {
  }

  // ACCESS FUNCTIONS
  Vertex_handle vertex(int i) const
  {
    CGAL_precondition( i >= 0 && i <= 3 );
    return V[i];
  }

  bool has_vertex(Vertex_handle v) const
  {
    return (V[0] == v) || (V[1] == v) || (V[2]== v) || (V[3]== v);
  }

  bool has_vertex(Vertex_handle v, int & i) const
  {
    if (v == V[0]) { i = 0; return true; }
    if (v == V[1]) { i = 1; return true; }
    if (v == V[2]) { i = 2; return true; }
    if (v == V[3]) { i = 3; return true; }
    return false;
  }

  int index(Vertex_handle v) const
  {
    if (v == V[0]) { return 0; }
    if (v == V[1]) { return 1; }
    if (v == V[2]) { return 2; }
    CGAL_assertion( v == V[3] );
    return 3;
  }

  Cell_handle neighbor(int i) const
  {
    CGAL_precondition( i >= 0 && i <= 3);
    return N[i];
  }

  bool has_neighbor(Cell_handle n) const
  {
    return (N[0] == n) || (N[1] == n) || (N[2] == n) || (N[3] == n);
  }

  bool has_neighbor(Cell_handle n, int & i) const
  {
    if(n == N[0]){ i = 0; return true; }
    if(n == N[1]){ i = 1; return true; }
    if(n == N[2]){ i = 2; return true; }
    if(n == N[3]){ i = 3; return true; }
    return false;
  }

  int index(Cell_handle n) const
  {
    if (n == N[0]) return 0;
    if (n == N[1]) return 1;
    if (n == N[2]) return 2;
    CGAL_assertion( n == N[3] );
    return 3;
  }

  // SETTING
  void set_neighbor(int i, Cell_handle n)
  {
    CGAL_precondition( i >= 0 && i <= 3);
    CGAL_precondition( this != n.operator->() );
    N[i] = n;
  }

  void set_neighbors()
  {
    N[0] = N[1] = N[2] = N[3] = Cell_handle();
  }

  void set_neighbors(Cell_handle n0, Cell_handle n1,
                     Cell_handle n2, Cell_handle n3)
  {
    CGAL_precondition( this != n0.operator->() );
    CGAL_precondition( this != n1.operator->() );
    CGAL_precondition( this != n2.operator->() );
    CGAL_precondition( this != n3.operator->() );
    N[0] = n0;
    N[1] = n1;
    N[2] = n2;
    N[3] = n3;
  }

  // CHECKING

  // the following trivial is_valid allows
  // the user of derived cell base classes
  // to add their own purpose checking
  bool is_valid(bool = false, int = 0) const
  { return true; }

  // For use by Compact_container.
  void * for_compact_container() const { return N[0].for_compact_container(); }
  void for_compact_container(void *p) { N[0].for_compact_container(p); }

  // TDS internal data access functions.
        TDS_data& tds_data()       { return _tds_data; }
  const TDS_data& tds_data() const { return _tds_data; }

  // We must override the functions that modify the vertices.
  // And if the point inside a vertex is modified, we fail,
  // but there's not much we can do for this now.
  void set_vertex(int i, Vertex_handle v)
  {
    CGAL_precondition( i >= 0 && i <= 3);
    V[i] = v;
  }

  void set_vertices()
  {
    V[0] = V[1] = V[2] = V[3] = Vertex_handle();
  }

  void set_vertices(Vertex_handle v0, Vertex_handle v1,
                    Vertex_handle v2, Vertex_handle v3)
  {
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    V[3] = v3;
  }

  // Returns the index of the cell of the input complex that contains the cell
  Subdomain_index subdomain_index() const { return subdomain_index_; }

  // Sets the index of the cell of the input complex that contains the cell
  void set_subdomain_index(const Subdomain_index& index)
  { subdomain_index_ = index; }

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

  /// Returns true if facet lies on a surface patch
  bool is_facet_on_surface(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return ( !( Surface_patch_index() == surface_index_table_[facet] ));
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
      Get_io_signature<Subdomain_index>()()
      + "+(" + Get_io_signature<Surface_patch_index>()() + ")[4]";
  }

  /// For the determinism of Compact_container iterators
  ///@{
  typedef Tag_true Has_timestamp;

  std::size_t time_stamp() const {
    return time_stamp_;
  }
  void set_time_stamp(const std::size_t& ts) {
    time_stamp_ = ts;
  }
  ///@}

private:
  /// Stores surface_index for each facet of the cell
  std::array<Surface_patch_index, 4> surface_index_table_ = {};

  std::array<Cell_handle, 4> N;
  std::array<Vertex_handle, 4> V;

  std::size_t time_stamp_;

  // The index of the cell of the input complex that contains me
  Subdomain_index subdomain_index_ = {};

  TDS_data _tds_data;

public:
  friend std::istream& operator>>(std::istream &is, Compact_simplicial_mesh_cell_3 &c)
  {
    Subdomain_index index;
    if(IO::is_ascii(is))
      is >> index;
    else
      read(is, index);

    if(is)
    {
      c.set_subdomain_index(index);
      for(int i = 0; i < 4; ++i)
      {
        Surface_patch_index i2;
        if(IO::is_ascii(is))
          is >> IO::iformat(i2);
        else
          read(is, i2);

        c.set_surface_patch_index(i, i2);
      }
    }

    return is;
  }

  friend
  std::ostream& operator<<(std::ostream &os, const Compact_simplicial_mesh_cell_3&c)
  {
    if(IO::is_ascii(os))
      os << c.subdomain_index();
    else
      write(os, c.subdomain_index());

    for(int i = 0; i < 4; ++i)
    {
      if(IO::is_ascii(os))
        os << ' ' << IO::oformat(c.surface_patch_index(i));
      else
        write(os, c.surface_patch_index(i));
    }

    return os;
  }
};

/*!
\ingroup PkgSMDS3Classes

The class `Compact_simplicial_mesh_cell_base_3`
is a model of the concept `SimplicialMeshCellBase_3`.
It is designed to serve as cell base class for.3D simplicial mesh data structures.
It stores and gives access to data about the complex the cell belongs to, such as the
subdomain it belongs to or surface it takes part to.
It is more compact in memory than `CGAL::Simplicial_mesh_cell_base_3`.

\tparam SubdomainIndex Type of indices for subdomains of the discretized geometric domain.
Must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible`
and `EqualityComparable`. The default constructed value must match the label
of the exterior of the domain (which contains at least the unbounded component).
It must match `MeshDomain_3::Subdomain_index` when used for mesh generation.

\tparam SurfacePatchIndex Type of indices for surface patches (boundaries and interfaces)
of the discretized geometric domain.
Must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible`
and `EqualityComparable`. The default constructed value must be the index value
assigned to a non surface facet.
It must match `MeshDomain_3::Surface_patch_index` when used for mesh generation.

\cgalModels{SimplicialMeshCellBase_3}

\sa `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveIndex>`
\sa \link Mesh_cell_base_3 `CGAL::Mesh_cell_base_3`\endlink
\sa `MeshDomain_3`
\sa `MeshDomainWithFeatures_3`
*/
template <typename SubdomainIndex,
          typename SurfacePatchIndex>
class Compact_simplicial_mesh_cell_base_3
{
public:
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type                              Triangulation_data_structure;
#else
  typedef internal::Dummy_tds_3                         Triangulation_data_structure;
#endif
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Cell_handle     Cell_handle;

public:
  template <typename TDS2>
  struct Rebind_TDS
  {
    typedef Compact_simplicial_mesh_cell_3<SubdomainIndex,
                                           SurfacePatchIndex,
                                           TDS2> Other;
  };
};

} // namespace CGAL

#endif // CGAL_COMPACT_SIMPLICIAL_MESH_CELL_BASE_3_H

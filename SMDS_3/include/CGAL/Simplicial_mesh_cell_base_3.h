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


#ifndef CGAL_SIMPLICIAL_MESH_CELL_BASE_3_H
#define CGAL_SIMPLICIAL_MESH_CELL_BASE_3_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/basic.h>
#include <CGAL/TDS_3/internal/Dummy_tds_3.h>
#include <CGAL/tags.h>
#include <CGAL/Has_timestamp.h>

#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/SMDS_3/io_signature.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <atomic>
#endif

namespace CGAL {

/*!
\ingroup PkgSMDS3Classes

The class `Simplicial_mesh_cell_base_3`
is a model of the concept `SimplicialMeshCellBase_3`.
It is designed to serve as cell base class for 3D simplicial mesh data structures.
It stores and gives access to data about the complex the cell belongs to, such as the
subdomain it belongs to or surface it takes part to.

\tparam Gt is the geometric traits class.
It must be a model of the concept `TriangulationTraits_3`

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

\tparam Cb is the cell base class from which `Simplicial_mesh_cell_base_3` derives.
It must be a model of the concept `TriangulationCellBase_3`.

\cgalModels{SimplicialMeshCellBase_3}

\sa `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveIndex>`
\sa \link Mesh_cell_base_3 `CGAL::Mesh_cell_base_3`\endlink
\sa `MeshDomain_3`
\sa `MeshDomainWithFeatures_3`
*/
template <typename Gt,
          typename SubdomainIndex,
          typename SurfacePatchIndex,
          typename Cb = Triangulation_cell_base_3<Gt> >
class Simplicial_mesh_cell_base_3
  : public Cb
{
public:
  using Vertex_handle = typename Cb::Vertex_handle;
  using Cell_handle = typename Cb::Cell_handle;

  using Geom_traits = Gt;

  // Index Type
  using Subdomain_index = SubdomainIndex;
  using Surface_patch_index = SurfacePatchIndex;

public:
  template <typename TDS2>
  struct Rebind_TDS
  {
    using Cb2 = typename Cb::template Rebind_TDS<TDS2>::Other;
    using Other = Simplicial_mesh_cell_base_3<Gt,
                                              SubdomainIndex,
                                              SurfacePatchIndex,
                                              Cb2>;
  };

public:
  Simplicial_mesh_cell_base_3()
    : time_stamp_(std::size_t(-1))
  {}

  Simplicial_mesh_cell_base_3(const Simplicial_mesh_cell_base_3& rhs)
    : Cb(static_cast<const Cb&>(rhs)),
      time_stamp_(rhs.time_stamp_),
      subdomain_index_(rhs.subdomain_index_)
  {
    for(int i=0; i <4; ++i)
      surface_index_table_[i] = rhs.surface_index_table_[i];
  }

  Simplicial_mesh_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                              Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3)
  { }

  Simplicial_mesh_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                              Vertex_handle v2, Vertex_handle v3,
                              Cell_handle n0, Cell_handle n1,
                              Cell_handle n2, Cell_handle n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3)
  { }

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

  std::size_t time_stamp_;

  // The index of the cell of the input complex that contains me
  Subdomain_index subdomain_index_ = {};

public:

  friend std::istream& operator>>(std::istream& is,
                                  Simplicial_mesh_cell_base_3& c)
  {
    Subdomain_index index;
    if(IO::is_ascii(is))
      is >> index;
    else
      read(is, index);

    if(is) {
      c.set_subdomain_index(index);
      for(int i = 0; i < 4; ++i)
      {
        Surface_patch_index i2;
        if(IO::is_ascii(is))
          is >> IO::iformat(i2);
        else
          {
            read(is, i2);
          }
        c.set_surface_patch_index(i, i2);
      }
    }
    return is;
  }

  friend
  std::ostream& operator<<(std::ostream& os,
                           const Simplicial_mesh_cell_base_3& c)
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

} // namespace CGAL

#endif // CGAL_SIMPLICIAL_MESH_CELL_BASE_3_H

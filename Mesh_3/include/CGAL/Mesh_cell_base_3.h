// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008 GeometryFactory, Sophia Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Stephane Tayeb


#ifndef CGAL_MESH_CELL_BASE_3_H
#define CGAL_MESH_CELL_BASE_3_H

#include <CGAL/license/Mesh_3.h>


// #define CGAL_DEPRECATED_HEADER "<CGAL/Mesh_cell_base_3.h>"
// #define CGAL_REPLACEMENT_HEADER "<CGAL/Compact_mesh_cell_base_3.h>"
// #include <CGAL/Installation/internal/deprecation_warning.h>

#include <CGAL/Mesh_3/config.h>

#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_with_weighted_circumcenter_3.h>
#include <CGAL/Mesh_3/Mesh_surface_cell_base_3.h>
#include <CGAL/SMDS_3/io_signature.h>
#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <atomic>
#endif

namespace CGAL {

// Sequential
template <typename Concurrency_tag>
class Mesh_cell_base_3_base
{
public:
#if defined(CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE) \
 || defined(CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE)

  // Erase counter (cf. Compact_container)
  unsigned int erase_counter() const
  {
    return this->m_erase_counter;
  }
  void set_erase_counter(unsigned int c)
  {
    this->m_erase_counter = c;
  }
  void increment_erase_counter()
  {
    ++this->m_erase_counter;
  }

private:
  typedef unsigned int              Erase_counter_type;
  Erase_counter_type                m_erase_counter;
#endif
};

#ifdef CGAL_LINKED_WITH_TBB
// Specialized version (parallel)
template <>
class Mesh_cell_base_3_base<Parallel_tag>
{
public:
  // Erase counter (cf. Compact_container)
  unsigned int erase_counter() const
  {
    return this->m_erase_counter;
  }
  void set_erase_counter(unsigned int c)
  {
    this->m_erase_counter = c;
  }
  void increment_erase_counter()
  {
    ++this->m_erase_counter;
  }

protected:
  typedef std::atomic<unsigned int> Erase_counter_type;
  Erase_counter_type                m_erase_counter;
};
#endif // CGAL_LINKED_WITH_TBB

// Class Mesh_cell_base_3
// Cell base class used in 3D meshing process.
// Adds information to Cb about the cell of the input complex containing it
/*!
\ingroup PkgMesh3MeshClasses
<!-- Meta-comment: this class cannot be deprecated by
Compact_mesh_cell_base_3, because the latter has a different API.
-- Laurent Rineau, 2013/10/16
\deprecated This class is deprecated since \cgal 4.3. Use
`CGAL::Compact_mesh_cell_base_3<GT,MD,Tds>` instead.
-->

The class `Mesh_cell_base_3<GT, MD, Cb>` is a model of the concept `MeshCellBase_3`.
It is designed to serve as cell base class for the 3D triangulation
used in the 3D mesh generation process.

\tparam GT is the geometric traits class.
It has to be a model of the concept `MeshTriangulationTraits_3`.

\tparam MD provides the types of indices used to identify
the faces of the input complex. It has to be a model
of the concept `MeshDomain_3`.

\tparam Cb is the cell base class. It has to be a model
of the concept `RegularTriangulationCellBaseWithWeightedCircumcenter_3` and defaults to
`Regular_triangulation_cell_base_with_weighted_circumcenter_3<GT>`.

\cgalModels{MeshCellBase_3}

\sa `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveIndex>`
\sa `CGAL::Compact_mesh_cell_base_3<GT, MD, Tds>`

*/
template<class GT,
         class MD,
         class Cb = CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3<
                      GT, CGAL::Regular_triangulation_cell_base_3<GT> > >
class Mesh_cell_base_3
#ifndef DOXYGEN_RUNNING
  : public Mesh_3::Mesh_surface_cell_base_3<GT, MD, Cb>,
    public Mesh_cell_base_3_base<typename Mesh_3::Mesh_surface_cell_base_3<GT, MD, Cb>::Tds::Concurrency_tag>
#endif
{
  typedef typename GT::FT FT;

public:
  // Base
  typedef Mesh_3::Mesh_surface_cell_base_3<GT, MD, Cb> Base;
  // Index Type
  typedef typename MD::Subdomain_index      Subdomain_index;
  typedef typename MD::Surface_patch_index  Surface_patch_index;

  // Backward compatibility
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index               Surface_index;
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX

  // Triangulation
  typedef typename Base::Tds Tds;
  typedef typename Tds::Vertex_handle Vertex_handle;
  typedef typename Tds::Cell_handle Cell_handle;

  // To get correct cell type in TDS
  template < class TDS3 >
  struct Rebind_TDS
  {
    typedef typename Cb::template Rebind_TDS<TDS3>::Other Cb3;
    typedef Mesh_cell_base_3 <GT, MD, Cb3> Other;
  };

  Mesh_cell_base_3()
    : Base()
    , subdomain_index_()
    , sliver_value_(FT(0.))
    , sliver_cache_validity_(false)
#ifdef CGAL_INTRUSIVE_LIST
    , next_intrusive_()
    , previous_intrusive_()
#endif
  {}

  Mesh_cell_base_3 (Vertex_handle v0,
                    Vertex_handle v1,
                    Vertex_handle v2,
                    Vertex_handle v3)
    : Base (v0, v1, v2, v3)
    , subdomain_index_()
    , sliver_value_(FT(0.))
    , sliver_cache_validity_(false)
#ifdef CGAL_INTRUSIVE_LIST
    , next_intrusive_()
    , previous_intrusive_()
#endif
  {}

  Mesh_cell_base_3 (Vertex_handle v0,
                    Vertex_handle v1,
                    Vertex_handle v2,
                    Vertex_handle v3,
                    Cell_handle n0,
                    Cell_handle n1,
                    Cell_handle n2,
                    Cell_handle n3)
    : Base (v0, v1, v2, v3, n0, n1, n2, n3)
    , subdomain_index_()
    , sliver_value_(FT(0.))
    , sliver_cache_validity_(false)
#ifdef CGAL_INTRUSIVE_LIST
    , next_intrusive_()
    , previous_intrusive_()
#endif
  {}

  // Default copy constructor and assignment operator are ok

  // Returns the index of the cell of the input complex that contains the cell
  Subdomain_index subdomain_index() const { return subdomain_index_; }

  // Sets the index of the cell of the input complex that contains the cell
  void set_subdomain_index(const Subdomain_index& index)
  { subdomain_index_ = index; }

  void set_sliver_value(const FT& value)
  {
    sliver_cache_validity_ = true;
    sliver_value_ = value;
  }
  const FT& sliver_value() const { return sliver_value_; }
  bool is_cache_valid() const { return sliver_cache_validity_; }
  void reset_cache_validity() const { sliver_cache_validity_ = false;  }

  /// \name I/O
  ///@{

  static
  std::string io_signature()
  {
    return
      Get_io_signature<Subdomain_index>()() + "+"
      + Get_io_signature<Base>()();
  }

  /// @}

#ifdef CGAL_INTRUSIVE_LIST
public:
  Cell_handle next_intrusive() const { return next_intrusive_; }
  void set_next_intrusive(Cell_handle c)
  {
    next_intrusive_ = c;
  }

  Cell_handle previous_intrusive() const { return previous_intrusive_; }
  void set_previous_intrusive(Cell_handle c)
  {
    previous_intrusive_ = c;
  }
#endif // CGAL_INTRUSIVE_LIST

  /// \name Determinism
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
  // The index of the cell of the input complex that contains me
  Subdomain_index subdomain_index_;

  FT sliver_value_;
  mutable bool sliver_cache_validity_;

#ifdef CGAL_INTRUSIVE_LIST
  Cell_handle next_intrusive_, previous_intrusive_;
#endif
  std::size_t time_stamp_;

};  // end class Mesh_cell_base_3

template < class GT, class MT, class Cb >
std::istream&
operator>>(std::istream &is,
           Mesh_cell_base_3<GT, MT, Cb> &c)
{
  typename Mesh_cell_base_3<GT, MT, Cb>::Subdomain_index index;
  if(IO::is_ascii(is))
    is >> index;
  else
    read(is, index);
  typedef typename Mesh_cell_base_3<GT, MT, Cb>::Base Cell_base;
  is >> static_cast<Cell_base&>(c);
  if(is) c.set_subdomain_index(index);
  return is;
}

template < class GT, class MT, class Cb >
std::ostream&
operator<<(std::ostream &os,
           const Mesh_cell_base_3<GT, MT, Cb> &c)
{
  if(IO::is_ascii(os))
     os << c.subdomain_index();
  else
    write(os, c.subdomain_index());
  typedef typename Mesh_cell_base_3<GT, MT, Cb>::Base Cell_base;
  return os << static_cast<const Cell_base&>(c);
}

}  // end namespace CGAL


#endif // CGAL_MESH_CELL_BASE_3_H

// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2011 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stéphane Tayeb, Andreas Fabri
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************


#ifndef CGAL_COMPACT_MESH_VERTEX_BASE_3_H
#define CGAL_COMPACT_MESH_VERTEX_BASE_3_H

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/internal/Mesh_3/indices_management.h>
#include <CGAL/Mesh_3/io_signature.h>
#include <CGAL/Has_timestamp.h>
#include <CGAL/tags.h>

namespace CGAL {

// Without erase counter
template <typename Concurrency_tag>
class Mesh_vertex_base_3_base
{
#if defined(CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE) \
 || defined(CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE)

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
  typedef unsigned int              Erase_counter_type;
  Erase_counter_type                m_erase_counter;
#endif
};

#ifdef CGAL_LINKED_WITH_TBB
// Specialized version (parallel)
template <>
class Mesh_vertex_base_3_base<Parallel_tag>
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
  typedef tbb::atomic<unsigned int> Erase_counter_type;
  Erase_counter_type                m_erase_counter;

};
#endif // CGAL_LINKED_WITH_TBB

// Class Mesh_vertex_base_3
// Vertex base class used in 3D meshing process.
// Adds information to Vb about the localization of the vertex in regards
// to the 3D input complex.
template<class GT,
         class Indices_tuple,
         class Index_,
         class Vb>
class Mesh_vertex_3
: public Vb,
  public Mesh_vertex_base_3_base<
    typename Vb::Triangulation_data_structure::Concurrency_tag>
{
public:
  typedef Vb Cmvb3_base;
  typedef typename Vb::Vertex_handle  Vertex_handle;

  // Types
  typedef Index_                      Index;
  typedef typename GT::FT             FT;

  // Constructor
  Mesh_vertex_3()
    : Vb()
    , number_of_incident_facets_(0)
    , number_of_components_(0)
    , index_()
    , meshing_info_(0)
    , dimension_(-1)
    , cache_validity(false)
#ifdef CGAL_INTRUSIVE_LIST
    , next_intrusive_()
    , previous_intrusive_()
#endif //CGAL_INTRUSIVE_LIST
  {}

  // Default copy constructor and assignment operator are ok

  // Returns the dimension of the lowest dimensional face of the input 3D
  // complex that contains the vertex
  int in_dimension() const {
    if(dimension_ < -1) return -2-dimension_;
    else return dimension_;
  }

  // Sets the dimension of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  void set_dimension(const int dimension) {
    CGAL_assertion(dimension < 4);
    dimension_ = short(dimension);
  }

  // Tells if the vertex is marked as a special protecting ball
  bool is_special() const { return dimension_ < -1; }

  // Marks or unmarks the vertex as a special protecting ball
  void set_special(bool special = true) {
    if(special != (dimension_ < -1) )
      dimension_ = short(-2-dimension_);
  }

  // Returns the index of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  Index index() const { return index_; }

  // Sets the index of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  void set_index(const Index& index) { index_ = index; }

  // Accessors to meshing_info private data
  const FT& meshing_info() const { return meshing_info_; }
  void set_meshing_info(const FT& value) { meshing_info_ = value; }

#ifdef CGAL_INTRUSIVE_LIST
  Vertex_handle next_intrusive() const { return next_intrusive_; }
  void set_next_intrusive(Vertex_handle v)
  {
    next_intrusive_ = v;
  }

  Vertex_handle previous_intrusive() const { return previous_intrusive_; }
  void set_previous_intrusive(Vertex_handle v)
  {
    previous_intrusive_ = v;
  }
#endif

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

  bool is_c2t3_cache_valid() const {
    return cache_validity;
  }

  void invalidate_c2t3_cache()
  {
    cache_validity = false;
  }

  void set_c2t3_cache(const std::size_t i, const std::size_t j)
  {
    number_of_incident_facets_ = i;
    number_of_components_ = j;
    cache_validity = true;
  }

  std::size_t cached_number_of_incident_facets() const
  {
    return number_of_incident_facets_;
  }

  std::size_t cached_number_of_components() const
  {
    return number_of_components_;
  }

  static
  std::string io_signature()
  {
    return
      Get_io_signature<Vb>()() + "+" +
      Get_io_signature<int>()() + "+" +
      Get_io_signature<Index>()();
  }
private:

  std::size_t number_of_incident_facets_;
  std::size_t number_of_components_; // number of components in the adjacency
  // graph of incident facets (in complex)


  // Index of the lowest dimensional face of the input 3D complex
  // that contains me
  Index index_;
  // Stores info needed by optimizers
  FT meshing_info_;

  // Dimension of the lowest dimensional face of the input 3D complex
  // that contains me. Negative values are a marker for special vertices.
  short dimension_;
  bool cache_validity;
#ifdef CGAL_INTRUSIVE_LIST
  Vertex_handle next_intrusive_;
  Vertex_handle previous_intrusive_;
#endif
  std::size_t time_stamp_;
public:

  friend std::istream& operator>>(std::istream &is, Mesh_vertex_3& v)
  {
    is >> static_cast<Cmvb3_base&>(v);
    int dimension;
    if(is_ascii(is)) {
      is >> dimension;

    } else {
      CGAL::read(is, dimension);
    }
    v.set_dimension(dimension);
    CGAL_assertion(v.in_dimension() >= -1);
    CGAL_assertion(v.in_dimension() < 4);
    Index index =
      Mesh_3::internal::Read_write_index<Indices_tuple,
                                         Index>()(is, v.in_dimension());
    v.set_index(index);
    return is;
  }

  friend std::ostream& operator<<(std::ostream &os, const Mesh_vertex_3& v)
  {
    os << static_cast<const Cmvb3_base&>(v);
    if(is_ascii(os)) {
      os << " " << v.in_dimension()
         << " ";
    } else {
      CGAL::write(os, v.in_dimension());
    }
    Mesh_3::internal::Read_write_index<Indices_tuple,
                                       Index>()(os,
                                                v.in_dimension(),
                                                v.index());
    return os;
  }
};  // end class Mesh_vertex_3

template<class GT,
         class MD,
         class Vb = Regular_triangulation_vertex_base_3<GT> >
struct Mesh_vertex_base_3 {
  using Triangulation_data_structure = internal::Dummy_tds_3;
  using Vertex_handle = typename Triangulation_data_structure::Vertex_handle;
  using Cell_handle = typename Triangulation_data_structure::Cell_handle;

  template < class TDS3 >
  struct Rebind_TDS {
    using Vb3 = typename Vb::template Rebind_TDS<TDS3>::Other;
    using Other = Mesh_vertex_3 <GT,
                                 Mesh_3::internal::Indices_tuple_t<MD>,
                                 typename MD::Index, Vb3>;
  };
};

template<class GT,
         class Indices_tuple,
         class Index,
         class Vb = Regular_triangulation_vertex_base_3<GT> >
struct Mesh_vertex_generator_3 {
  using Triangulation_data_structure = internal::Dummy_tds_3;
  using Vertex_handle = typename Triangulation_data_structure::Vertex_handle;
  using Cell_handle = typename Triangulation_data_structure::Cell_handle;

  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS3>::Other Vb3;
    typedef Mesh_vertex_3 <GT, Indices_tuple, Index, Vb3> Other;
  };
};

}  // end namespace CGAL



#endif // CGAL_COMPACT_MESH_VERTEX_BASE_3_H

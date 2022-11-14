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
// Author(s)     : Stéphane Tayeb, Andreas Fabri, Jane Tournois
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************


#ifndef CGAL_SIMPLICIAL_MESH_VERTEX_BASE_3_H
#define CGAL_SIMPLICIAL_MESH_VERTEX_BASE_3_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/Triangulation_vertex_base_3.h>

#include <CGAL/SMDS_3/internal/indices_management.h>
#include <CGAL/SMDS_3/io_signature.h>
#include <CGAL/Has_timestamp.h>

#include <boost/variant.hpp>
#include <tuple>

namespace CGAL {

//Class Simplicial_mesh_vertex_3
//Adds information to Vb about the localization of the vertex in regards
//to the 3D input complex.
//\cgalModels `SimplicialMeshVertexBase_3`
template<class GT,
         typename Subdomain_index,
         typename Surface_patch_index,
         typename Curve_index,
         typename Corner_index,
         class Vb>
class Simplicial_mesh_vertex_3
: public Vb
{
private :
  using Cmvb3_base = Vb;

public:
  using Vertex_handle = typename Vb::Vertex_handle;

  // Types
  using Index = boost::variant<Subdomain_index, Surface_patch_index, Curve_index, Corner_index>;
  using FT = typename GT::FT;

  // Constructor
  Simplicial_mesh_vertex_3()
    : Vb()
    , number_of_incident_facets_(0)
    , number_of_components_(0)
    , index_()
    , dimension_(-1)
    , cache_validity(false)
    , time_stamp_(std::size_t(-1))
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

  // Returns the index of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  Index index() const { return index_; }

  // Sets the index of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  void set_index(const Index& index) { index_ = index; }

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

  // Dimension of the lowest dimensional face of the input 3D complex
  // that contains me. Negative values are a marker for special vertices.
  short dimension_;
  bool cache_validity;
  std::size_t time_stamp_;

public:

  friend std::istream& operator>>(std::istream &is, Simplicial_mesh_vertex_3& v)
  {
    is >> static_cast<Cmvb3_base&>(v);
    int dimension;
    if(IO::is_ascii(is)) {
      is >> dimension;

    } else {
      CGAL::read(is, dimension);
    }
    v.set_dimension(dimension);
    CGAL_assertion(v.in_dimension() >= -1);
    CGAL_assertion(v.in_dimension() < 4);

    using Indices_tuple = std::tuple<Subdomain_index,
                                     Surface_patch_index,
                                     Curve_index,
                                     Corner_index>;
    Index index =
      Mesh_3::internal::Read_write_index<Indices_tuple,
                                         Index>()(is, v.in_dimension());
    v.set_index(index);
    return is;
  }

  friend std::ostream& operator<<(std::ostream &os, const Simplicial_mesh_vertex_3& v)
  {
    os << static_cast<const Cmvb3_base&>(v);
    if(IO::is_ascii(os)) {
      os << " " << v.in_dimension()
         << " ";
    } else {
      CGAL::write(os, v.in_dimension());
    }
    using Indices_tuple = std::tuple<Subdomain_index,
                                     Surface_patch_index,
                                     Curve_index,
                                     Corner_index>;
    Mesh_3::internal::Read_write_index<Indices_tuple,
                                       Index>()(os,
                                                v.in_dimension(),
                                                v.index());
    return os;
  }
};  // end class Simplicial_mesh_vertex_3


/*!
\ingroup PkgSMDS3Classes

The class `Simplicial_mesh_vertex_base_3` is a model of the concept
`SimplicialMeshVertexBase_3`.
It is designed to serve as vertex base class for 3D simplicial mesh data structures.
It stores and gives access to data about the complex the vertex belongs to, such as the
index of the subcomplex it belongs to.

\tparam Gt is the geometric traits class.
It must be a model of the concept `TriangulationTraits_3`

\tparam Subdomain_index Type of indices for subdomains of the discretized geometric domain.
Must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible`
and `EqualityComparable`. The default constructed value must match the label
of the exterior of the domain (which contains at least the unbounded component).
 It must match the `Subdomain_index` of the model
  of the `MeshDomain_3` concept when used for mesh generation.

\tparam Surface_patch_index Type of indices for surface patches (boundaries and interfaces)
of the discretized geometric domain.
Must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible`
and `EqualityComparable`. The default constructed value must be the index value
assigned to a non surface facet.
 It must match the `Surface_patch_index` of the model
  of the `MeshDomain_3` concept when used for mesh generation.

\tparam Curve_index Type of indices for curves (i.e. \f$ 1\f$-dimensional features)
of the discretized geometric domain.
Must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible` and
`LessThanComparable`. The default constructed value must be the value for an edge which
does not approximate a 1-dimensional feature of the geometric domain.
 It must match the `Curve_index` types of the model
  of the `MeshDomainWithFeatures_3` concept when used for mesh generation.

\tparam  Corner_index Type of indices for corners (i.e.\f$ 0\f$--dimensional features)
of the discretized geometric domain.
It must be a model of `CopyConstructible`, `Assignable`, `DefaultConstructible` and
`LessThanComparable`.
 It must match the `Corner_index` of the model
  of the `MeshDomainWithFeatures_3` concept when used for mesh generation.

\tparam Vb is the vertex base class. It has to be a model
of the concept `TriangulationVertexBase_3` and defaults to
`Triangulation_vertex_base_3<Gt>`.

\cgalModels `SimplicialMeshVertexBase_3`

\sa `CGAL::Mesh_complex_3_in_triangulation_3`
\sa \link Mesh_vertex_base_3 `CGAL::Mesh_vertex_base_3`\endlink
\sa `MeshDomain_3`
\sa `MeshDomainWithFeatures_3`
*/
template<class Gt,
         typename Subdomain_index,
         typename Surface_patch_index,
         typename Curve_index,
         typename Corner_index,
         class Vb = Triangulation_vertex_base_3<Gt> >
struct Simplicial_mesh_vertex_base_3
{
  using Triangulation_data_structure = internal::Dummy_tds_3;
  using Vertex_handle = typename Triangulation_data_structure::Vertex_handle;
  using Cell_handle = typename Triangulation_data_structure::Cell_handle;

  template < class TDS3 >
  struct Rebind_TDS {
    using Vb3 = typename Vb::template Rebind_TDS<TDS3>::Other;
    using Other = Simplicial_mesh_vertex_3 <Gt,
      Subdomain_index,
      Surface_patch_index,
      Curve_index,
      Corner_index,
      Vb3>;
  };
};


}  // end namespace CGAL



#endif // CGAL_SIMPLICIAL_MESH_VERTEX_BASE_3_H

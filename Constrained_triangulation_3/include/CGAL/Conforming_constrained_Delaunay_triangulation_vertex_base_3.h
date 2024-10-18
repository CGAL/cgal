// Copyright (c) 2019-2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/Base_with_time_stamp.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_data_3.h>
#include <CGAL/SMDS_3/io_signature.h>

namespace CGAL {

/**
 * @ingroup PkgCT_3Classes
 * @brief The Conforming_constrained_Delaunay_triangulation_vertex_base_3 class is a vertex base class for the
 *        conforming constrained Delaunay triangulation in 3D.
 *
 * This class is derived from the `Triangulation_vertex_base_3` class and provides additional functionality
 * required by `make_conforming_constrained_Delaunay_triangulation_3()`.
 *
 * @tparam Traits The geometric traits class, model of `DelaunayTriangulationTraits_3`.
 *         It must be the same as the geometric traits class of the triangulation.
 * @tparam Vertex_base The base class for the vertex. It must be a model of `TriangulationVertexBase_3`.
 *
 * @cgalModels{ConformingConstrainedDelaunayTriangulationVertexBase_3}
 *
 * \sa `CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3`
 */
template < typename Traits, typename Vertex_base = Triangulation_vertex_base_3<Traits> >
class Conforming_constrained_Delaunay_triangulation_vertex_base_3 : public Base_with_time_stamp<Vertex_base>
{
  Conforming_constrained_Delaunay_triangulation_vertex_data_3 ccdt_3_data_;
  bool cache_validity_ = false;
  CDT_3_face_index index_ = 0;
  int dim_ = -1;
  std::size_t number_of_incident_facets_ = 0;
  std::size_t number_of_components_ = 0;

public:
  // To get correct vertex type in TDS
  template <class TDS3> struct Rebind_TDS
  {
    using Vb3 = typename Vertex_base::template Rebind_TDS<TDS3>::Other;
    using Other = Conforming_constrained_Delaunay_triangulation_vertex_base_3<Traits, Vb3>;
  };

  // constructors, inherited from the base class
  using Base = Base_with_time_stamp<Vertex_base>;
  using Base::Base;

  // model of SimplicialMeshVertexBase_3
  using Index = CDT_3_face_index;
  int in_dimension() const { return dim_; }
  void set_dimension(int d) { dim_ = d; }
  Index index() const { return index_; }
  void set_index(Index i) { index_ = i; }
  bool is_c2t3_cache_valid() const { return cache_validity_; }
  void invalidate_c2t3_cache() { cache_validity_ = false; }
  void set_c2t3_cache(std::size_t i, std::size_t j)
  {
    number_of_incident_facets_ = i;
    number_of_components_ = j;
    cache_validity_ = true;
  }
  std::size_t cached_number_of_incident_facets() const { return number_of_incident_facets_; }
  std::size_t cached_number_of_components() const { return number_of_components_; }

  void sync() {
    switch(ccdt_3_data().vertex_type()) {
      case CDT_3_vertex_type::FREE:
        set_dimension(3);
        set_index(0);
        break;
      case CDT_3_vertex_type::CORNER:
        set_dimension(0);
        set_index(0);
        break;
      case CDT_3_vertex_type::STEINER_ON_EDGE:
        set_dimension(1);
        set_index(0);
        break;
      case CDT_3_vertex_type::STEINER_IN_FACE:
        set_dimension(2);
        set_index(ccdt_3_data().face_index());
        break;
      default:
        CGAL_error();
        break;
    }
  }

  // model of ConformingConstrainedDelaunayTriangulationVertexBase_3
  Conforming_constrained_Delaunay_triangulation_vertex_data_3& ccdt_3_data() { return ccdt_3_data_; }
  const Conforming_constrained_Delaunay_triangulation_vertex_data_3& ccdt_3_data() const { return ccdt_3_data_; }

  static std::string io_signature() {
    return Get_io_signature<Vertex_base>()();
  }
};

} // namespace CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_BASE_3_H

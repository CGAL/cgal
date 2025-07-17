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

#include <CGAL/Triangulation_simplex_base_with_time_stamp.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_data_3.h>
#include <CGAL/SMDS_3/io_signature.h>

namespace CGAL {

/**
 * @ingroup PkgConstrainedTriangulation3Classes
 * @brief Vertex base class for the 3D conforming constrained Delaunay triangulation.
 *
 * This class is derived from its parameter template `VertexBase` and provides additional functionality
 * required by `Conforming_constrained_Delaunay_triangulation_3`.
 *
 * @tparam Traits The geometric traits class, model of `ConformingConstrainedDelaunayTriangulationTraits_3`.
 *         It must be the same as the geometric traits class of the triangulation.
 * @tparam VertexBase The base class for the vertex. It must be a model of `TriangulationVertexBase_3`.
 *
 * @cgalModels{ConformingConstrainedDelaunayTriangulationVertexBase_3}
 *
 * \sa `CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3`
 */
template < typename Traits, typename VertexBase = Triangulation_vertex_base_3<Traits> >
class Conforming_constrained_Delaunay_triangulation_vertex_base_3
  : public Triangulation_simplex_base_with_time_stamp<VertexBase>
{
  Conforming_constrained_Delaunay_triangulation_vertex_data_3 ccdt_3_data_;

public:
  // To get correct vertex type in TDS
  template <class TDS3> struct Rebind_TDS
  {
    using Vb3 = typename VertexBase::template Rebind_TDS<TDS3>::Other;
    using Other = Conforming_constrained_Delaunay_triangulation_vertex_base_3<Traits, Vb3>;
  };

  // constructors, inherited from the base class
  using Base = Triangulation_simplex_base_with_time_stamp<VertexBase>;
  using Base::Base;

  // model of ConformingConstrainedDelaunayTriangulationVertexBase_3
  Conforming_constrained_Delaunay_triangulation_vertex_data_3& ccdt_3_data() { return ccdt_3_data_; }
  const Conforming_constrained_Delaunay_triangulation_vertex_data_3& ccdt_3_data() const { return ccdt_3_data_; }

  static std::string io_signature() {
    return Get_io_signature<VertexBase>()();
  }
};

} // namespace CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_BASE_3_H

// Copyright (c) 2019-2024  GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_vertex_data_3.h>
#include <CGAL/SMDS_3/io_signature.h>

namespace CGAL {

/**
 * @brief The Constrained_Delaunay_triangulation_vertex_base_3 class is a vertex base class for the
 *        Constrained Delaunay Triangulation in 3D.
 *
 * This class is derived from the `Triangulation_vertex_base_3` class and provides additional functionality
 * required by `Constrained_Delaunay_triangulation_3`.
 *
 * @tparam Gt The geometric traits class, model of `DelaunayTriangulationTraits_3`.
 *         It must be the same as the geometric traits class of the triangulation.
 * @tparam Vb The base class for the vertex. It must be a model of `TriangulationVertexBase_3`.
 *
 * @cgalModels{ConstrainedDelaunayTriangulationVertexBase_3}
 *
 * \sa `CGAL::Constrained_Delaunay_triangulation_cell_base_3`
 */
template < typename Gt, typename Vb = Triangulation_vertex_base_3<Gt> >
class Constrained_Delaunay_triangulation_vertex_base_3 : public Base_with_time_stamp<Vb>
{
  Constrained_Delaunay_triangulation_vertex_data_3 cdt_3_data_;
public:
  // To get correct vertex type in TDS
  template <class TDS3> struct Rebind_TDS
  {
    using Vb3 = typename Vb::template Rebind_TDS<TDS3>::Other;
    using Other = Constrained_Delaunay_triangulation_vertex_base_3<Gt, Vb3>;
  };

  using Base = Base_with_time_stamp<Vb>;
  using Base::Base;

  Constrained_Delaunay_triangulation_vertex_data_3& cdt_3_data() {
    return cdt_3_data_;
  }

  const Constrained_Delaunay_triangulation_vertex_data_3& cdt_3_data() const {
    return cdt_3_data_;
  }

  static std::string io_signature() {
    return Get_io_signature<Vb>()();
  }
};

} // namespace CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_BASE_3_H

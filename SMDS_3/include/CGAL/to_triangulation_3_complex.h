// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_TO_TRIANGULATION_3_COMPLEX_H
#define CGAL_TO_TRIANGULATION_3_COMPLEX_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Triangulation_3.h>

namespace CGAL
{

/*!
 * \ingroup PkgSMDS3Functions
 * converts a `Mesh_complex_3_in_triangulation_3`
 * to a `Mesh_complex_3_in_triangulation_3` with a `Triangulation_3` as
 * underlying triangulation.
 *
 * @tparam C3t3 model of `MeshComplex_3InTriangulation_3`
 * @tparam CornerIndex is the type of the corner indices
 * @tparam CurveIndex is the type of the curve indices
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @todo
 */
template <typename C3t3,
          typename NamedParameters = parameters::Default_named_parameters>
CGAL::Mesh_complex_3_in_triangulation_3<
  CGAL::Triangulation_3<typename C3t3::Triangulation::Geom_traits,
                        typename C3t3::Triangulation::Triangulation_data_structure>>
to_triangulation_3_complex(
  C3t3 c3t3,
  const NamedParameters& np = parameters::default_values())
{
  using GT = typename C3t3::Triangulation::Geom_traits;
  using TDS = typename C3t3::Triangulation::Triangulation_data_structure;
  using Corner_index = typename C3t3::Corner_index;
  using Curve_index = typename C3t3::Curve_index;

  using T3 = CGAL::Triangulation_3<GT, TDS>;
  using C3t3_new = CGAL::Mesh_complex_3_in_triangulation_3<T3, Corner_index, Curve_index>;

  T3 tr;
  tr.swap(c3t3.triangulation());
  C3t3_new c3t3_new;
  c3t3_new.triangulation().swap(tr);

  // todo : port metadata from input c3t3 to the new one

  return c3t3_new;
}

} // end namespace CGAL

#endif // CGAL_TO_TRIANGULATION_3_COMPLEX_H

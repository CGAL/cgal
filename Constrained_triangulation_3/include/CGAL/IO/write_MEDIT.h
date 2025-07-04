// Copyright (c) 2019-2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jane Tournois

#include "CGAL/unordered_flat_map.h"
#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/IO/File_medit.h>

#include <stack>
#include <ostream>

namespace CGAL
{
namespace IO
{
/*!
 * @ingroup PkgCDT3IOFunctions
 * @brief outputs a conforming constrained Delaunay triangulation to
 * the MEDIT (`.mesh`) file format.
 *        See \cgalCite{frey:inria-00069921} for a comprehensive description of this
 *        file format.
 * @param os the output stream
 * @param ccdt the conforming constrained Delaunay triangulation to be written
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamNBegin{with_plc_face_id}
 *    \cgalParamDescription{a Boolean activating the numbering of PLC face identifiers in the output.
 *                          If `ccdt` was constructed with the `plc_face_id` property map given as a named parameter,
 *                          and this parameter is set to `true`,
 *                          the output will contain the corresponding patch identifier for each facet of the triangulation.
 *                          If this parameter is set to `false`, the output will not contain any patch identifier.
 *                          If `ccdt` was not constructed with the `plc_face_id` property map, and this parameter is
 *                          set to `true`, the output will contain a patch identifier for each facet of the triangulation.}
 *    \cgalParamType{Boolean}
 *    \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd

 * \see \ref IOStreamMedit
 */
template <typename Traits,
          typename Tr,
          typename NamedParameters = parameters::Default_named_parameters>
void write_MEDIT(std::ostream& os,
                 const Conforming_constrained_Delaunay_triangulation_3<Traits, Tr>& ccdt,
                 const NamedParameters& np = parameters::default_values())
{
  const auto& tr = ccdt.triangulation();

  using Tr_ = typename cpp20::remove_cvref_t<decltype(tr)>;

  using Vertex_handle = typename Tr_::Vertex_handle;
  using Facet = typename Tr_::Facet;
  using Cell_handle = typename Tr_::Cell_handle;
  CGAL::unordered_flat_map<Cell_handle, bool> cells_in_domain;
  for(Cell_handle c : tr.all_cell_handles())
    cells_in_domain[c] = true;

  std::stack<Cell_handle> stack;
  stack.push(tr.infinite_cell());
  while(!stack.empty())
  {
    auto ch = stack.top();
    stack.pop();
    cells_in_domain[ch] = false;
    for(int i = 0; i < 4; ++i)
    {
      if(ccdt.is_facet_constrained(ch, i))
        continue;
      auto n = ch->neighbor(i);
      if(cells_in_domain[n])
        stack.push(n);
    }
  }

  using parameters::choose_parameter;
  using parameters::get_parameter;
  const bool has_plc_face_id
    = choose_parameter(get_parameter(np, internal_np::with_plc_face_id), false);

  auto plc_patch_map = boost::make_function_property_map<Facet>([&](const Facet& f)
    { return has_plc_face_id ? f.first->ccdt_3_data().face_constraint_index(f.second) + 1 : 1; });

  return SMDS_3::output_to_medit(os,
                                 tr,
                                 tr.finite_vertex_handles(),
                                 ccdt.constrained_facets(),
                                 tr.finite_cell_handles(),
                                 boost::make_function_property_map<Vertex_handle>([](Vertex_handle) { return 0; }),
                                 plc_patch_map,
                                 boost::make_function_property_map<Cell_handle>([&](Cell_handle ch) {
                                   return cells_in_domain[ch] ? 1 : 0;
                                 }));
}

}// end namespace IO
}// end namespace CGAL

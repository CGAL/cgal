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

#include "CGAL/type_traits.h"
#include "CGAL/unordered_flat_map.h"
#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/IO/File_medit.h>
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
 *
 * \see \ref IOStreamMedit
 */
template <typename Traits, typename Tr_of_default>
void write_MEDIT(std::ostream& os,
                 const Conforming_constrained_Delaunay_triangulation_3<Traits, Tr_of_default>& ccdt)
{
  const auto& tr = ccdt.triangulation();

  using Tr = typename cpp20::remove_cvref_t<decltype(tr)>;

  using Vertex_handle = typename Tr::Vertex_handle;
  using Facet = typename Tr::Facet;
  using Cell_handle = typename Tr::Cell_handle;
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

  return SMDS_3::output_to_medit(os,
                                 tr,
                                 tr.finite_vertex_handles(),
                                 ccdt.constrained_facets(),
                                 tr.finite_cell_handles(),
                                 boost::make_function_property_map<Vertex_handle>([](Vertex_handle) { return 0; }),
                                 boost::make_function_property_map<Facet>([](const Facet& f) {
                                   return f.first->ccdt_3_data().face_constraint_index(f.second) + 1;
                                 }),
                                 boost::make_function_property_map<Cell_handle>([&](Cell_handle ch) {
                                   return cells_in_domain[ch] ? 1 : 0;
                                 }));
}

}// end namespace IO
}// end namespace CGAL

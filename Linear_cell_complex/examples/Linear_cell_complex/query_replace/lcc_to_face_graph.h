// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef CGAL_LCC_TO_FACE_GRAPH_H
#define CGAL_LCC_TO_FACE_GRAPH_H

#include <unordered_map>
#include <vector>
#include <CGAL/config.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/property_map.h>
#include <boost/unordered_map.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/iterator/function_output_iterator.hpp>

namespace CGAL {

/** Build a face graph from an LCC. Only 3-free faces are used in order to
 *  construct in tm only the surface of lcc. */
template <typename LCC, typename TargetMesh>
void lcc_to_face_graph(const LCC& lcc, TargetMesh& tm)
{
  typedef typename LCC::Dart_const_handle DH;
  typedef typename LCC::Vertex_attribute_const_handle VH;

  typedef typename boost::graph_traits<TargetMesh>::vertex_descriptor tm_vertex_descriptor;
  
  std::unordered_map<VH, tm_vertex_descriptor> vertices_map;
  std::vector<tm_vertex_descriptor> vertices;

  for (auto it=lcc.vertex_attributes().begin();
       it!=lcc.vertex_attributes().end(); ++it)
  {
    tm_vertex_descriptor vd=CGAL::add_vertex(tm);
    tm.point(vd)=it->point();
    vertices_map[it]=vd;
  }

  auto treated=lcc.get_new_mark();
  for (auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  {
    if (lcc.template is_free<3>(it) && !lcc.is_marked(it, treated))
    {
      vertices.clear();
      DH cur=it;
      do
      {
        vertices.push_back(vertices_map[lcc.vertex_attribute(cur)]);
        lcc.mark(cur, treated);
        cur=lcc.next(cur);
      }
      while(cur!=it);

      if (!CGAL::Euler::can_add_face(vertices, tm))
      { std::cerr<<"[ERROR] a face cannot be added !!"<<std::endl; }

      CGAL::Euler::add_face(vertices, tm);
    }
    else
    { lcc.mark(it, treated); }
  }

  assert(lcc.is_whole_map_marked(treated));
  lcc.free_mark(treated);
}

} // namespace CGAL

#endif //  CGAL_LCC_TO_FACE_GRAPH_H
////////////////////////////////////////////////////////////////////////////////

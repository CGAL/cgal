// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_BGL_PARTITION_H
#define CGAL_BGL_PARTITION_H

/**
* \ingroup PkgBGL
* \file CGAL/boost/graph/partition.h
* Convenience header file including the headers for all the partitioning-related
* free functions of this package.
*/

#include <CGAL/boost/graph/METIS/partition_graph.h>
#include <CGAL/boost/graph/METIS/partition_dual_graph.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace CGAL {

namespace internal {

// Note that to use the function below with Polyhedron_3, you need to enhance
// the Polyhedron_3 to use items, that is, use:
// typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> PM;

// \ingroup PkgBGLPartition
//
// Output each part of a partition as a single mesh.
//
// \param tm a triangle mesh
// \param nparts the number of parts
// \param fpmap the property map with the partition indices
// \param filename_base Partitions will be output in `.off` files named
//                      `{filename_base}_[0...nparts].off`
//
// \tparam TriangleMesh must be a model of a `FaceListGraph`, `HalfedgeListGraph`, and \bgllink{VertexListGraph}.
// \tparam FacePartitionIDPmap is a model of `ReadablePropertyMap`
//           with `boost::graph_traits<TriangleMesh>::%face_descriptor`
//           as key type and `boost::graph_traits<Graph>::%faces_size_type` as value type.
template<typename TriangleMesh, typename FacePartitionIDPmap>
void output_partition(const TriangleMesh& tm,
                      const idx_t nparts,
                      const FacePartitionIDPmap fpmap,
                      const std::string filename_base)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  typedef CGAL::Face_filtered_graph<TriangleMesh>         Filtered_graph;

  for(int i=0; i<nparts; ++i)
  {
    std::ostringstream filename;
    filename << filename_base << "_" << i << ".off" << std::ends;
    std::ofstream out(filename.str().c_str());

    Filtered_graph m_part(tm, i, fpmap);
    if(!m_part.is_selection_valid())
    {
      std::cerr << "Warning: cannot extract subdomain #" << i << " because it is "
                << "not a manifold mesh" << std::endl;
      continue;
    }

    TriangleMesh out_mesh;
    CGAL::copy_face_graph(m_part, out_mesh);

    out << out_mesh;
  }
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_BGL_PARTITION_H

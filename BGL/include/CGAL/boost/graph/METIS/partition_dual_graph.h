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

#ifndef CGAL_BGL_PARTITION_DUAL_GRAPH_H
#define CGAL_BGL_PARTITION_DUAL_GRAPH_H

#include <CGAL/boost/graph/METIS/partition_graph.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/assertions.h>

#include <metis.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/tuple/tuple.hpp>

namespace CGAL {

namespace METIS {

template<typename TriangleMesh, typename METIS_options, typename NamedParameters>
void partition_dual_graph(const TriangleMesh& tm, int nparts,
                          METIS_options options, // options array
                          const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));
  CGAL_precondition_msg(nparts > 1, ("Partitioning requires a number of parts > 1"));

  using boost::choose_param;
  using boost::get_param;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_iterator       face_iterator;

  // vertex index map
  typedef typename CGAL::Polygon_mesh_processing::GetVertexIndexMap<TriangleMesh, NamedParameters>::type Indices;
  Indices indices = choose_param(get_param(np, internal_np::vertex_index),
                                 get_const_property_map(boost::vertex_index, tm));

  idx_t nn = static_cast<idx_t>(num_vertices(tm));
  idx_t ne = static_cast<idx_t>(num_faces(tm));
  idx_t d = 3; // number of nodes per element
  idx_t* eptr = new idx_t[ne + 1];
  idx_t* eind = new idx_t[d * ne];

  // fill the adjacency info
  face_iterator fit, fe;
  boost::tie(fit, fe) = faces(tm);
  for(int i=0, j=0; fit!=fe; ++fit, ++i)
  {
    eptr[i] = j;

    halfedge_descriptor h = halfedge(*fit, tm), done = h;
    do
    {
      vertex_descriptor v = target(h, tm);
      CGAL_assertion(j < d * ne);
      eind[j++] = static_cast<idx_t>(get(indices, v));
      h = next(h, tm);
    } while (h != done);

    CGAL_assertion(i < ne);
    eptr[i + 1] = j;
  }

  // a dual edge between elements exists if they share 'nparts' vertices
  idx_t ncommon = 2;

  // either the edgecut or the total communication volume of the dual graph’s partitioning
  idx_t objval;

  // partition info for the nodes
  idx_t* npart = (idx_t*) calloc(nn, sizeof(idx_t));
  CGAL_assertion(npart != NULL);

  // partition info for the elements
  idx_t* epart = (idx_t*) calloc(ne, sizeof(idx_t));
  CGAL_assertion(epart != NULL);

  // do not support Fortran-style arrays
  CGAL_assertion((*options)[METIS_OPTION_NUMBERING] == -1 || // default initialization is '-1'
                 (*options)[METIS_OPTION_NUMBERING] == 0);

  CGAL_assertion_code(int ret =)
    METIS_PartMeshDual(&ne, &nn, eptr, eind,
                       NULL /* elements weights*/, NULL /*elements sizes*/,
                       &ncommon, &nparts,
                       NULL /* partitions weights */,
                       *options,
                       &objval, epart, npart);

  CGAL_assertion(ret == METIS_OK);

  Output_vertex_partition_ids vo;
  Output_face_partition_ids fo;
  vo(tm, indices, npart, get_param(np, internal_np::vertex_partition_id));
  fo(tm, epart, get_param(np, internal_np::face_partition_id));
}

template<typename TriangleMesh, typename NamedParameters>
void partition_dual_graph(const TriangleMesh& tm, int nparts,
                          const boost::param_not_found, // no METIS options were passed
                          const NamedParameters& np)
{
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  return partition_dual_graph(tm, nparts, &options, np);
}

/// \ingroup PkgBGLPartition
///
/// Computes a partition of the input triangular mesh into `nparts` parts,
/// based on the mesh's dual graph. The resulting partition is stored in the vertex and/or face
/// property maps that are passed as parameters using \ref bgl_namedparameters "Named Parameters".
///
/// Property map for `CGAL::vertex_index_t` should be either available
/// as an internal property map to `tm` or provided as \ref bgl_namedparameters "Named Parameters".
///
/// \param tm a triangle mesh
/// \param nparts the number of parts in the final partition
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \tparam TriangleMesh is a model of the `FaceListGraph` concept.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_index_map}
///     is a property map containing the index of each vertex of `tm` intialized from `0` to `num_vertices(tm)-1`.
///   \cgalParamEnd
///   \cgalParamBegin{METIS_options}
///     is a parameter used in to pass options to the METIS mesh
///     partitioner. The many options of METIS are not described here. Instead, users
///     should refer to the <a href="http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf">documentation</a>
///     of METIS directly.
///   \cgalParamEnd
///   \cgalParamBegin{vertex_partition_id_map}
///     is a property map that contains (after the function has been run)
///     the ID of the subpart for each vertex of `tm`.
///   \cgalParamEnd
///   \cgalParamBegin{face_partition_id_map}
///     is a property map that contains (after the function has been run)
///     the ID of the subpart for each face of `tm`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \pre `tm` is a pure triangular surface mesh: there are no edges
///       without at least one incident face
template<typename TriangleMesh, typename NamedParameters>
void partition_dual_graph(const TriangleMesh& tm, int nparts, const NamedParameters& np)
{
  using boost::get_param;

  return partition_dual_graph(tm, nparts, get_param(np, internal_np::METIS_options), np);
}

template<typename TriangleMesh>
void partition_dual_graph(const TriangleMesh& tm, const int nparts)
{
  return partition_dual_graph(tm, nparts, CGAL::parameters::all_default());
}

} // end namespace METIS

} // end namespace CGAL

#endif // CGAL_BGL_PARTITION_DUAL_GRAPH_H

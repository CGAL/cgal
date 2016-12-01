// Copyright (c) 2016  GeometryFactory (France).
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
//
//
// Author(s)     :

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONE_HELPER_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONE_HELPER_H

#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/boost/graph/properties.h>

#include <boost/foreach.hpp>
#include <boost/graph/graph_traits.hpp>

#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace CGAL {

namespace Surface_mesh_parameterization {

enum Cone_type
{
  Unique_cone = 0,
  Duplicated_cone
};

namespace internal {

/// Read the cones from the input file.
template<typename Polygon_mesh>
Error_code read_cones(const Polygon_mesh& pm, const char* filename,
  std::vector<typename boost::graph_traits<Polygon_mesh>::vertex_descriptor>& cone_vds_in_pm)
{
  typedef typename boost::graph_traits<Polygon_mesh>::vertex_descriptor PM_vertex_descriptor;

  std::ifstream in(filename);
  std::string vertices_line;
  std::getline(in, vertices_line); // read the first line of the file
  std::istringstream iss(vertices_line);
  std::vector<int> cones;
  cones.reserve(4);
  int cone_index;
  while(iss >> cone_index) {
    cones.push_back(cone_index);
  }

  if(cones.size() < 3 || cones.size() > 4) {
    std::cerr << "Error: Not enough or too many input cones" << std::endl;
    return ERROR_WRONG_PARAMETER;
  }

  std::cout << "Cones: ";
  for(std::size_t i=0; i<cones.size(); ++i)
    std::cout << cones[i] << " ";
  std::cout << std::endl;

  // Locate the cones in the underlying mesh 'pm'
  CGAL_assertion(cone_vds_in_pm.empty());
  cone_vds_in_pm.resize(cones.size());

  for(std::size_t i=0; i<cones.size(); ++i) {
    int counter = 0;
    BOOST_FOREACH(PM_vertex_descriptor vd, vertices(pm)) {
      if(counter == cones[i]) {
        cone_vds_in_pm[i] = vd;
        break;
      }
      ++counter;
    }
    CGAL_postcondition(cone_vds_in_pm[i] != PM_vertex_descriptor());
  }

  return OK;
}

/// Locate the cones on the seam mesh (find the corresponding seam mesh
/// vertex_descriptor) and mark them with a tag that indicates whether it is a
/// simple cone or a duplicated cone.
///
/// The cones are ordered: the first and last cones are the extremetities of the seam.
///
/// \tparam Mesh is a seam mesh
/// \tparam BaseMesh is the underlying mesh of `Mesh`
/// \tparam ConeMap a map vertex_descriptor --> Cone_type
template<typename Mesh,
         typename BaseMesh,
         typename Cones_in_pmesh_vector,
         typename ConeMap>
void locate_cones(const Mesh& mesh,
                  const Cones_in_pmesh_vector& cone_vds_in_sm,
                  ConeMap& cones)
{
  CGAL_precondition(cones.empty());

  typedef typename boost::graph_traits<BaseMesh>::vertex_descriptor   BM_vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor       vertex_descriptor;

  // property map to go from BM_vertex_descriptor to Point_3
  typedef typename Kernel_traits<BaseMesh>::PPM        PM_PPM;
  const PM_PPM pm_ppmap = get(boost::vertex_point, mesh.mesh());

  // property map to go from vertex_descriptor to Point_3
  typedef typename Kernel_traits<Mesh>::PPM            PPM;
  const PPM ppmap = get(boost::vertex_point, mesh);

  // the cones in the underlying mesh
  std::size_t cvdss = cone_vds_in_sm.size();
  for(std::size_t i=0; i<cvdss; ++i) {
    BM_vertex_descriptor smvd = cone_vds_in_sm[i];
    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
      if(get(ppmap, vd) == get(pm_ppmap, smvd)) {
        Cone_type ct = (i == 0 || i == cvdss-1) ? Unique_cone : Duplicated_cone;
        cones.insert(std::make_pair(vd, ct));
      }
    }
  }

  CGAL_postcondition((cone_vds_in_sm.size() == 3 && cones.size() == 4) ||
                     (cone_vds_in_sm.size() == 4 && cones.size() == 6));

  std::cout << cone_vds_in_sm.size() << " cones in sm" << std::endl;
  std::cout << cones.size() << " cones in mesh" << std::endl;
}

/// Same as above, but the cones are NOT ordered and we thus use seam mesh
/// information to determine which cones are Duplicate_cones and which cones
/// are unique
///
/// \tparam Mesh is a seam mesh
/// \tparam BaseMesh is the type of the underlying mesh in the seam mesh
/// \tparam ConeMap a map vertex_descriptor --> Cone_type
template<typename Mesh,
         typename BaseMesh,
         typename Cones_in_pmesh_set,
         typename ConeMap>
void locate_unordered_cones(const Mesh& mesh,
                            const Cones_in_pmesh_set& cone_vds_in_sm,
                            ConeMap& cones)
{
  CGAL_precondition(cones.empty());

  typedef typename boost::graph_traits<BaseMesh>::vertex_descriptor     BM_vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor       halfedge_descriptor;

  // find a vertex on the seam
  vertex_descriptor vertex_on_seam;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
    if(mesh.has_on_seam(vd)) {
      vertex_on_seam = vd;
      break;
    }
  }

  CGAL_assertion(vertex_on_seam != vertex_descriptor());

  // property map to go from BM_vertex_descriptor to Point_3
  typedef typename Kernel_traits<BaseMesh>::PPM        PM_PPM;
  const PM_PPM pm_ppmap = get(boost::vertex_point, mesh.mesh());

  // property map to go from vertex_descriptor to Point_3
  typedef typename Kernel_traits<Mesh>::PPM            PPM;
  const PPM ppmap = get(boost::vertex_point, mesh);

  // walk on the seam and mark if we encounter a cone
  vertex_descriptor end = vertex_on_seam;
  do {
    BOOST_FOREACH(BM_vertex_descriptor smvd, cone_vds_in_sm) {
      if(get(ppmap, vertex_on_seam) == get(pm_ppmap, smvd)) { // the seam mesh vertex is a cone

        // we have encountered a cone. Must check if the cone is a Unique_cone
        // or a Duplicated_cone.

        // a check is to look at the sources of two halfedges with the same direction
        // on either side of the seam. If the sources are the same, it's a Unique_cone;
        // if they differ, it's a duplicated_cone.

        halfedge_descriptor hd = halfedge(vertex_on_seam, mesh);
        halfedge_descriptor other_hd = opposite(hd, mesh);

        // little trick to go from border halfedge on one side of the seam to
        // interior halfedge on the other side of the seam
        CGAL_assertion(other_hd.seam);
        other_hd.seam = false;

        Cone_type ct = (target(hd, mesh) == source(other_hd, mesh)) ? Unique_cone
                                                                    : Duplicated_cone;

        std::cout << "new cone with type: " << ct << std::endl;

        cones.insert(std::make_pair(vertex_on_seam, ct));
      }
    }

    // move to the next vertex_descriptor on the seam
    vertex_on_seam = source(halfedge(vertex_on_seam, mesh), mesh);
    CGAL_assertion(mesh.has_on_seam(vertex_on_seam));

  } while(vertex_on_seam != end);

  CGAL_postcondition((cone_vds_in_sm.size() == 3 && cones.size() == 4) ||
                     (cone_vds_in_sm.size() == 4 && cones.size() == 6));

  std::cout << cone_vds_in_sm.size() << " cones in sm" << std::endl;
  std::cout << cones.size() << " cones in mesh" << std::endl;
}

} // namespace internal

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONE_HELPER_H

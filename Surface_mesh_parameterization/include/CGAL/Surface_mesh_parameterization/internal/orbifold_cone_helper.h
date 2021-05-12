// Copyright (c) 2016  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONE_HELPER_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONE_HELPER_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/orbifold_enums.h>

#include <CGAL/boost/graph/properties.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_set.hpp>

#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

/// \file orbifold_cone_helper.h

namespace CGAL {

namespace Surface_mesh_parameterization {

namespace internal {

// Orbifold type functions
template<typename Point_container>
Point_container get_cones_parameterized_coordinates(const Orbifold_type orb_type)
{
  typedef typename Point_container::value_type                Point;

  Point_container tcoords;
  if(orb_type == Square) {
    tcoords.push_back(Point(-1, -1));
    tcoords.push_back(Point(1, 1));
  } else if(orb_type == Parallelogram) {
    tcoords.push_back(Point(0, -0.5));
    tcoords.push_back(Point(-1, -0.5));
    tcoords.push_back(Point(0, 0.5));
  } else { // if(orb_type == Diamond || orb_type == Triangle)
    tcoords.push_back(Point(-1, 1));
    tcoords.push_back(Point(-1, -1));
  }

  return tcoords;
}

template<typename NT_container>
NT_container get_angles_at_cones(const Orbifold_type orb_type)
{
  // Note that angles are minus as we go around the seam border in a counterclockwise manner
  NT_container angs;
  if(orb_type == Square) {
    angs.push_back(4.);
    angs.push_back(4.);
  } else if(orb_type == Diamond) {
    angs.push_back(3.);
    angs.push_back(3.);
  } else if(orb_type == Triangle) {
    angs.push_back(6.);
    angs.push_back(2.);
  } else { // if(orb_type == Parallelogram)
    angs.push_back(2);
    angs.push_back(1);
    angs.push_back(2);
  }

  return angs;
}

template<typename Cone_container>
bool are_cones_unique(const Cone_container& cones)
{
  std::size_t n_of_cones = cones.size();
  if(n_of_cones == 0) {
    std::cerr << "Warning: Check cone uniquess with no cones...?" << std::endl;
    return true;
  }

  typedef typename Cone_container::value_type   Cone;
  boost::unordered_set<Cone> unique_cones;

  unique_cones.insert(cones.begin(), cones.end());

  return (n_of_cones == unique_cones.size());
}

// Locate the cone tagged 'First_unique_cone' and its index in the seam mesh.
template<typename vertex_descriptor,
         typename ConeMap,
         typename VertexIndexMap>
void find_start_cone(const ConeMap& cmap,
                     VertexIndexMap vimap,
                     vertex_descriptor& cone,
                     int& cone_index)
{
  CGAL_precondition(!cmap.empty());

  typename ConeMap::const_iterator cmit = cmap.begin(), cend = cmap.end();
  for(; cmit!=cend; ++cmit) {
    if(cmit->second != First_unique_cone)
      continue;

    cone = cmit->first;
    cone_index = get(vimap, cone);

    return;
  }

  std::cerr << "Error: did not find first cone" << std::endl;
  CGAL_postcondition(false);
}

// Locate the cone tagged 'First_unique_cone' in the seam mesh.
template<typename vertex_descriptor,
         typename ConeMap>
void find_start_cone(const ConeMap& cmap,
                     vertex_descriptor& cone)
{
  CGAL_precondition(!cmap.empty());

  typename ConeMap::const_iterator cmit = cmap.begin(), cend = cmap.end();
  for(; cmit!=cend; ++cmit) {
    if(cmit->second != First_unique_cone)
      continue;

    cone = cmit->first;
    return;
  }

  std::cerr << "Error: did not find first cone" << std::endl;
  CGAL_postcondition(false);
}

/// Check the validity of the input cones in the `Seam_mesh` mesh.
template<typename SeamMesh,
         typename ConeInputBidirectionalIterator,
         typename Cones_in_Seam_mesh_map>
bool check_cone_validity(const SeamMesh& mesh,
                         ConeInputBidirectionalIterator first, ConeInputBidirectionalIterator beyond,
                         const Cones_in_Seam_mesh_map& cones)
{
  typedef typename boost::graph_traits<SeamMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<SeamMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<typename SeamMesh::TriangleMesh>::vertex_descriptor TM_vertex_descriptor;

  typename std::iterator_traits<
    ConeInputBidirectionalIterator>::difference_type number_of_cones_in_tm = std::distance(first, beyond);

  // check cone numbers
  if((number_of_cones_in_tm == 3 && cones.size() != 4) ||
     (number_of_cones_in_tm == 4 && cones.size() != 6)) {
    std::cerr << "Invalid cone placement: " << std::endl;
    std::cerr << number_of_cones_in_tm << " cones in the base mesh" << std::endl;
    std::cerr << cones.size() << " cones in the seam mesh" << std::endl;
    return false;
  }

  // check cone types
  bool found_first_unique_cone = false, found_second_unique_cone = false;
  int duplicated_cone_counter = 0;
  typename Cones_in_Seam_mesh_map::const_iterator it = cones.begin(),
                                                  end = cones.end();
  for(; it!=end; ++it) {
    if(it->second == First_unique_cone) {
      if(found_first_unique_cone) {
        std::cerr << "Error: More than one 'First_unique_cone'" << std::endl;
        return false;
      }
      found_first_unique_cone = true;
    }
    else if(it->second == Second_unique_cone) {
      if(found_second_unique_cone) {
        std::cerr << "Error: More than one 'Second_unique_cone'" << std::endl;
        return false;
      }
      found_second_unique_cone = true;
    }
    else {
      if(it->second != Duplicated_cone) {
        std::cerr << "Error: Unknow cone type: " << it->second << std::endl;
        return false;
      }
      ++duplicated_cone_counter;
    }
  }

  if(!found_first_unique_cone || !found_second_unique_cone) {
    std::cerr << "Error: Could not find all unique cones" << std::endl;
    return false;
  }

  if((number_of_cones_in_tm == 3 && duplicated_cone_counter != 2) ||
     (number_of_cones_in_tm == 4 && duplicated_cone_counter != 4)) {
    std::cerr << "Error: Wrong number of duplicated cones" << std::endl;
    return false;
  }

  // check seams
  vertex_descriptor first_cone;
  find_start_cone(cones, first_cone);

  halfedge_descriptor hd = halfedge(first_cone, mesh);
  CGAL_precondition(mesh.has_on_seam(hd));
  halfedge_descriptor bhd = opposite(hd, mesh);
  CGAL_precondition(is_border(bhd, mesh));

  // count how many times vertices on a seam appear
  boost::unordered_map<TM_vertex_descriptor, int> seam_vertices_counter;

  for(halfedge_descriptor hdaf : halfedges_around_face(bhd, mesh)) {
    CGAL_precondition(mesh.has_on_seam(hdaf));
    TM_vertex_descriptor tm_vds = source(hdaf, mesh.mesh());
    TM_vertex_descriptor tm_vdt = target(hdaf, mesh.mesh());

    // insert vds
    std::pair<typename boost::unordered_map<TM_vertex_descriptor, int>::iterator,
              bool> is_insert_successful =
                      seam_vertices_counter.insert(std::make_pair(tm_vds, 1));
    if(!is_insert_successful.second)
      ++(is_insert_successful.first->second);

    // insert vdt
    is_insert_successful = seam_vertices_counter.insert(std::make_pair(tm_vdt, 1));
    if(!is_insert_successful.second)
      ++(is_insert_successful.first->second);
  }

  // check that the seam forms one connected component
  if(seam_vertices_counter.size() != (mesh.number_of_seam_edges() + 1)) {
    std::cerr << "Seam is not a connected component" << std::endl;
    return false;
  }

  // check for self intersections in the seam
  typename boost::unordered_map<TM_vertex_descriptor, int>::iterator sit = seam_vertices_counter.begin(),
                                                                     send = seam_vertices_counter.end();
  for(; sit!=send; ++sit) {
    if(sit->second != 2 && sit->second != 4) {
      std::cerr << sit->second << std::endl;
      std::cerr << "Seam intersect itself (or something bad like that...)" << std::endl;
      return false;
    }
  }

  // check that unique vertices are actually unique and duplicated cones are
  // actually duplicated
  it = cones.begin();
  for(; it!=end; ++it) {
    vertex_descriptor vd = it->first;
    TM_vertex_descriptor cvd = target(vd.hd, mesh.mesh());
    sit = seam_vertices_counter.find(cvd);
    if(sit == seam_vertices_counter.end()) {
      std::cerr << "Cone not on the seam" << std::endl;
      return false;
    }

    if(it->second == First_unique_cone || it->second == Second_unique_cone) {
      if(sit->second != 2) {
        std::cerr << "First or second cone not at the beginning/end of the seam" << std::endl;
        return false;
      }
    } else if(it->second == Duplicated_cone) {
      if(sit->second != 4) {
        std::cerr << "Duplicate cone not in the interior of the seam" << std::endl;
        return false;
      }
    }
  }

  return true;
}

} // namespace internal

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONE_HELPER_H

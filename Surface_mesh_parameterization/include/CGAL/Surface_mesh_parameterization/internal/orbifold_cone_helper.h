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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONE_HELPER_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONE_HELPER_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/boost/graph/properties.h>

#include <boost/foreach.hpp>
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

/// \ingroup PkgSurfaceParameterizationEnums
///
/// The types of cones used in Orbifold Tutte parameterization.
///
/// `Unique_cones` are found at the beginning and the end of the seam. All other
/// cones are duplicated in the sense that when the seam is `opened`, the vertex
/// is duplicated at two different positions.
enum Cone_type
{
  First_unique_cone = 0,
  Second_unique_cone,
  Duplicated_cone
};

/// \ingroup PkgSurfaceParameterizationEnums
///
/// The four Orbifold types available in the Orbifold Tutte parameterization.
enum Orbifold_type
{
  Square = 0,
  Diamond,
  Triangle,
  Parallelogram
};

/// Get message corresponding to an error code
/// \param orb_type The integer value in the enum
/// \return         The string describing the Orbifold type
const char* get_orbifold_type(int orb_type)
{
  // Messages corresponding to Error_code list above. Must be kept in sync!
  static const char* type[Parallelogram+1] = {
    "Square",
    "Diamond",
    "Triangle",
    "Parallelogram"
  };

  if(orb_type > Parallelogram || orb_type < 0)
    return "Unknown orbifold type";
  else
    return type[orb_type];
}

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

namespace internal {

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
         typename Cones_in_Seam_mesh_map,
         typename Cones_in_Base_mesh_container>
bool check_input_validity(const SeamMesh& mesh,
                          const Cones_in_Seam_mesh_map& cones,
                          const Cones_in_Base_mesh_container& cone_tm_vds)
{
  typedef typename boost::graph_traits<SeamMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<SeamMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<typename SeamMesh::TriangleMesh>::vertex_descriptor TM_vertex_descriptor;

  // check cone numbers
  if((cone_tm_vds.size() == 3 && cones.size() != 4) ||
     (cone_tm_vds.size() == 4 && cones.size() != 6)) {
    std::cerr << "Error: Problem in number of cones: " << std::endl;
    std::cerr << cone_tm_vds.size() << " cones in the base mesh" << std::endl;
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

  if((cone_tm_vds.size() == 3 && duplicated_cone_counter != 2) ||
     (cone_tm_vds.size() == 4 && duplicated_cone_counter != 4)) {
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

  BOOST_FOREACH(halfedge_descriptor hdaf, halfedges_around_face(bhd, mesh)) {
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

/// Read the cones from an input stream.
template<typename TriangleMesh>
Error_code read_cones(const TriangleMesh& pm, std::ifstream& in,
  std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cone_vds_in_tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor TM_vertex_descriptor;

  std::vector<int> cones;
  cones.reserve(4);
  int cone_index;
  while(in >> cone_index) {
    cones.push_back(cone_index);
  }

  std::cout << "Input cones: ";
  for(std::size_t i=0; i<cones.size(); ++i)
    std::cout << cones[i] << " ";
  std::cout << std::endl;

  if(cones.size() < 3 || cones.size() > 4) {
    std::cerr << "Error: Not enough or too many input cones" << std::endl;
    return ERROR_WRONG_PARAMETER;
  }

  if(!are_cones_unique(cones)) {
    std::cerr << "Error: The input cones are not unique" << std::endl;
    return ERROR_WRONG_PARAMETER;
  }

  // Locate the cones in the underlying mesh 'pm'
  CGAL_assertion(cone_vds_in_tm.empty());
  cone_vds_in_tm.resize(cones.size());

  for(std::size_t i=0; i<cones.size(); ++i) {
    int counter = 0;
    BOOST_FOREACH(TM_vertex_descriptor vd, vertices(pm)) {
      if(counter == cones[i]) {
        cone_vds_in_tm[i] = vd;
        break;
      }
      ++counter;
    }
    CGAL_postcondition(cone_vds_in_tm[i] != TM_vertex_descriptor());
  }

  return OK;
}

/// Read the cones from an input file.
template<typename TriangleMesh>
Error_code read_cones(const TriangleMesh& pm, const char* filename,
  std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cone_vds_in_tm)
{
  std::ifstream in(filename);
  return read_cones(pm, in, cone_vds_in_tm);
}


/// Locate the cones on the seam mesh (find the corresponding seam mesh
/// vertex_descriptor) and mark them with a tag that indicates whether it is a
/// simple cone or a duplicated cone.
///
/// The cones are ordered: the first and last cones are the extremetities of the seam.
///
/// \tparam SeamMesh is a seam mesh
/// \tparam ConeMap a map vertex_descriptor --> Cone_type
template<typename SeamMesh,
         typename Cones_in_pmesh_vector,
         typename ConeMap>
bool locate_cones(const SeamMesh& mesh,
                  const Cones_in_pmesh_vector& cone_tm_vds,
                  ConeMap& cones)
{
  CGAL_precondition(cones.empty());

  typedef typename SeamMesh::TriangleMesh                                  TriangleMesh;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    TM_vertex_descriptor;
  typedef typename boost::graph_traits<SeamMesh>::vertex_descriptor        vertex_descriptor;

  // property map to go from TM_vertex_descriptor to Point_3
  typedef typename Kernel_traits<TriangleMesh>::PPM                        PM_PPM;
  const PM_PPM pm_ppmap = get(boost::vertex_point, mesh.mesh());

  // property map to go from vertex_descriptor to Point_3
  typedef typename Kernel_traits<SeamMesh>::PPM                            PPM;
  const PPM ppmap = get(boost::vertex_point, mesh);

  // the cones in the underlying mesh
  std::size_t cvdss = cone_tm_vds.size();

  for(std::size_t i=0; i<cvdss; ++i) {
    TM_vertex_descriptor smvd = cone_tm_vds[i];
    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
      if(get(ppmap, vd) == get(pm_ppmap, smvd)) {
        Cone_type ct;
        if(i == 0)
          ct = First_unique_cone;
        else if(i == (cvdss - 1))
          ct = Second_unique_cone;
        else
          ct = Duplicated_cone;

        cones.insert(std::make_pair(vd, ct));
      }
    }
  }

  return check_input_validity(mesh, cones, cone_tm_vds);
}

/// Same as above, but the cones are NOT ordered and we thus use seam mesh
/// information to determine which cones are Duplicate_cones and which cones
/// are unique
///
/// \tparam SeamMesh is a seam mesh
/// \tparam Cones_in_pmesh_set is a set of cones (vertex_descriptor of SeamMesh)
/// \tparam ConeMap a map vertex_descriptor of SeamMesh  --> Cone_type
template<typename SeamMesh,
         typename Cones_in_pmesh_set,
         typename ConeMap>
bool locate_unordered_cones(const SeamMesh& mesh,
                            const Cones_in_pmesh_set& cone_tm_vds,
                            ConeMap& cones)
{
  CGAL_precondition(cones.empty());
  CGAL_precondition(cone_tm_vds.size() == 3 || cone_tm_vds.size() == 4);

  typedef typename SeamMesh::TriangleMesh                                  TriangleMesh;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    TM_vertex_descriptor;
  typedef typename boost::graph_traits<SeamMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<SeamMesh>::halfedge_descriptor      halfedge_descriptor;

  // find a vertex on the seam
  vertex_descriptor vertex_on_seam;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
    if(mesh.has_on_seam(vd)) {
      vertex_on_seam = vd;
      break;
    }
  }

  CGAL_assertion(vertex_on_seam != vertex_descriptor());

  // property map to go from TM_vertex_descriptor to Point_3
  typedef typename Kernel_traits<TriangleMesh>::PPM          PM_PPM;
  const PM_PPM pm_ppmap = get(boost::vertex_point, mesh.mesh());

  // property map to go from vertex_descriptor to Point_3
  typedef typename Kernel_traits<SeamMesh>::PPM              PPM;
  const PPM ppmap = get(boost::vertex_point, mesh);

  bool first_cone_met = false;

  // walk on the seam and mark if we encounter a cone
  vertex_descriptor end = vertex_on_seam;
  do {
    BOOST_FOREACH(TM_vertex_descriptor smvd, cone_tm_vds) {
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

        Cone_type ct;
        if(target(hd, mesh) == source(other_hd, mesh)) {
          if(first_cone_met)
            ct = Second_unique_cone;
          else {
            ct = First_unique_cone;
            first_cone_met = true;
          }
        } else {
          ct = Duplicated_cone;
        }

        cones.insert(std::make_pair(vertex_on_seam, ct));
      }
    }

    // move to the next vertex_descriptor on the seam
    vertex_on_seam = source(halfedge(vertex_on_seam, mesh), mesh);
    CGAL_assertion(mesh.has_on_seam(vertex_on_seam));

  } while(vertex_on_seam != end);

  return check_input_validity(mesh, cones, cone_tm_vds);
}

} // namespace internal

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONE_HELPER_H

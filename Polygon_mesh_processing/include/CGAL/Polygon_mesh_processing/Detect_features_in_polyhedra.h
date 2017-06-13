// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent Rineau, Stephane Tayeb
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYHEDRA_H
#define CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYHEDRA_H

#include <CGAL/license/Polygon_mesh_processing.h>


#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Polygon_mesh_processing/Detect_features_in_polyhedra_fwd.h>
#include <CGAL/Compare_handles_with_or_without_timestamps.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Mesh_3/properties.h>
#include <set>

namespace CGAL {
namespace Polygon_mesh_processing {

  template <typename Polyhedron, typename FT, typename Patch_id>
void detect_features(Polyhedron& p,
                     FT angle_in_deg,
                     typename boost::property_map<Polyhedron, CGAL::face_patch_id_t<Patch_id> >::type)
{
  Detect_features_in_polyhedra<Polyhedron, Patch_id> go;
  // AF todo: Add overload for the next three functions so that we use the pid_map
  //          Add a default for pid_map
  go.detect_sharp_edges(p,angle_in_deg);
  go.detect_surface_patches(p);
  go.detect_vertices_incident_patches(p);
}

  
  template <typename Polyhedron_, typename Patch_id_>
class Detect_features_in_polyhedra
{
public:

  typedef Polyhedron_ Polyhedron;
  typedef Patch_id_ Patch_id;

  typedef typename boost::property_traits<typename boost::property_map<Polyhedron,
  vertex_point_t>::type>::value_type Point_3;

  typedef typename Kernel_traits<Point_3>::Kernel Geom_traits;

  typedef typename Geom_traits::Vector_3    Vector_3;
  typedef typename Geom_traits::FT          FT;

  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor  vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor     face_descriptor;
  
  typedef typename boost::property_map<Polyhedron, face_patch_id_t<Patch_id> >::type PatchId_pmap;

  typedef CGAL::Compare_handles_with_or_without_timestamps Compare_handles;
  
  typedef std::set<face_descriptor> face_descriptor_set;
  typedef std::set<halfedge_descriptor> He_handle_set;
  
public:
  Detect_features_in_polyhedra()
    : current_surface_index_(1)
  {}
  
  void detect_sharp_edges(Polyhedron& polyhedron,
                          FT angle_in_deg = FT(60)) const;
  void detect_surface_patches(Polyhedron& polyhedron);
  void detect_vertices_incident_patches(Polyhedron& p);

  int maximal_surface_patch_index() const {
    return current_surface_index_ - 1;
  }
  
private:
  Vector_3 facet_normal(const Polyhedron& polyhedron, const face_descriptor& f) const;
  bool is_sharp(Polyhedron& polyhedron, const halfedge_descriptor& he, FT cos_angle) const;
  void flood(Polyhedron& polyhedron,
             face_descriptor f, const Patch_id& id,
             face_descriptor_set& unsorted_faces) const;

  template <typename Int>
  Int generate_patch_id(Int, int);

  template <typename Int>
  std::pair<Int, Int> generate_patch_id(std::pair<Int, Int>, int);
  
private:
  // Stores the current surface index (usefull to detect different patches
  // on different polyhedra)
  int current_surface_index_;
};


template <typename P_, typename I_>
void
Detect_features_in_polyhedra<P_, I_>::
detect_sharp_edges(Polyhedron& polyhedron, FT angle_in_deg) const
{
  // Initialize vertices
  typename boost::property_map<Polyhedron,vertex_num_feature_edges_t>::type vnfe
    = get(vertex_num_feature_edges_t(),polyhedron);

  typename boost::property_map<Polyhedron,halfedge_is_feature_t>::type hif
    = get(halfedge_is_feature_t(),polyhedron);

  BOOST_FOREACH(typename boost::graph_traits<P_>::vertex_descriptor vd, vertices(polyhedron))
  {
    put(vnfe,vd, 0);
  }
  
  FT cos_angle ( std::cos(CGAL::to_double(angle_in_deg) * CGAL_PI / 180.) );
  
  // Detect sharp edges
  BOOST_FOREACH(typename boost::graph_traits<P_>::edge_descriptor ed, edges(polyhedron))
  {
    typename boost::graph_traits<P_>::halfedge_descriptor he = halfedge(ed,polyhedron);
    if(is_border(he,polyhedron) || angle_in_deg == FT() ||
       (angle_in_deg != FT(180) && is_sharp(polyhedron,he,cos_angle))
       )
    {
      put(hif, he, true);
      put(hif, opposite(he,polyhedron), true);
      
      put(vnfe, target(he,polyhedron), get(vnfe, target(he,polyhedron))+1);
      put(vnfe, source(he,polyhedron), get(vnfe, source(he,polyhedron))+1);
    }
  }
}


template <typename P_, typename I_>
template <typename Int>
Int
Detect_features_in_polyhedra<P_, I_>::
generate_patch_id(Int, int i)
{
  return Int(i);
}

template <typename P_, typename I_>
template <typename Int>
std::pair<Int, Int>
Detect_features_in_polyhedra<P_, I_>::
generate_patch_id(std::pair<Int, Int>, int i)
{
  return std::pair<Int, Int>(i, 0);
}

template <typename P_, typename I_>
void
Detect_features_in_polyhedra<P_, I_>::
detect_surface_patches(Polyhedron& polyhedron)
{
  PatchId_pmap pid_map = get(face_patch_id_t<Patch_id>(),polyhedron);
  // Initialize unsorted_faces
  face_descriptor_set unsorted_faces;
  BOOST_FOREACH(typename boost::graph_traits<Polyhedron>::face_descriptor fd, faces(polyhedron))
  {
    unsorted_faces.insert(fd);
  }
  
  // Flood
  while ( ! unsorted_faces.empty() )
  {
    face_descriptor f = *(unsorted_faces.begin());
    unsorted_faces.erase(unsorted_faces.begin());
    
    const Patch_id patch_id = generate_patch_id(Patch_id(),
                                                current_surface_index_);
    put(pid_map, f, patch_id);
    flood(polyhedron, f,patch_id,unsorted_faces);
    ++current_surface_index_;
  }
}


template <typename P_, typename I_>
void
Detect_features_in_polyhedra<P_, I_>::
detect_vertices_incident_patches(Polyhedron& polyhedron)
{
  PatchId_pmap pid_map = get(face_patch_id_t<Patch_id>(),polyhedron);
  typename boost::property_map<Polyhedron,halfedge_is_feature_t>::type hif
      = get(halfedge_is_feature,polyhedron);
  typedef typename boost::property_map<Polyhedron,vertex_incident_patches_t<Patch_id> >::type VIP_map;
  VIP_map vip = get(vertex_incident_patches_t<Patch_id>(),polyhedron);

  BOOST_FOREACH(vertex_descriptor vit,vertices(polyhedron))
  {
    // Look only at feature vertices
    if( ! get(hif, halfedge(vit, polyhedron)) ){ continue; }
    
    // Loop on incident facets of vit
    std::set<Patch_id> set;
    BOOST_FOREACH(halfedge_descriptor he, halfedges_around_target(vit,polyhedron))
    {
      if( ! is_border(he,polyhedron) )
      {
        set.insert(get(pid_map,face(he,polyhedron)));
      }
      else if( ! is_border(opposite(he,polyhedron),polyhedron) )
      {
        set.insert(get(pid_map, face(opposite(he,polyhedron),polyhedron)));
      }
    }
    put(vip, vit, set);
  }
}
  
// -----------------------------------
// Private methods
// -----------------------------------
template <typename P_, typename I_>
typename Detect_features_in_polyhedra<P_, I_>::Vector_3
Detect_features_in_polyhedra<P_, I_>::
facet_normal(const Polyhedron& polyhedron, const face_descriptor& f) const
{
  return Polygon_mesh_processing::compute_face_normal(f,polyhedron);
}


template <typename P_, typename I_>
bool
Detect_features_in_polyhedra<P_, I_>::
is_sharp(Polyhedron& polyhedron, const halfedge_descriptor& he, FT cos_angle) const
{
  if(is_border(edge(he,polyhedron),polyhedron)){
    return false;
  }
  face_descriptor f1 = face(he,polyhedron);
  face_descriptor f2 = face(opposite(he,polyhedron),polyhedron);

  const Vector_3& n1 = facet_normal(polyhedron, f1);
  const Vector_3& n2 = facet_normal(polyhedron, f2);
  
  if ( n1 * n2 <= cos_angle )
    return true;
  else
    return false;
}
  
 
template <typename P_, typename I_>
void
Detect_features_in_polyhedra<P_, I_>::
flood(Polyhedron& polyhedron,
      face_descriptor f,
      const Patch_id &patch_id,
      face_descriptor_set& unsorted_faces) const
{
  PatchId_pmap pid_map = get(face_patch_id_t<Patch_id>(),polyhedron);
  typename boost::property_map<Polyhedron,halfedge_is_feature_t>::type hif
      = get(halfedge_is_feature,polyhedron);
  // Initialize he_to_explore with halfedges of the starting facet
  He_handle_set he_to_explore;
  BOOST_FOREACH(halfedge_descriptor hd,
                halfedges_around_face(halfedge(f,polyhedron), polyhedron))
  {
    he_to_explore.insert(opposite(hd,polyhedron));
  }
  
  // While there is something to explore
  while ( ! he_to_explore.empty() )
  {
    // Get next halfedge to explore
    halfedge_descriptor he = *(he_to_explore.begin());
    he_to_explore.erase(he_to_explore.begin());
    
    // If we don't go through a border of the patch
    if ( ! get(hif, he) && ! is_border(he,polyhedron) )
    {
      face_descriptor explored_facet = face(he,polyhedron);
      
      // Mark facet and delete it from unsorted
      put(pid_map, explored_facet, patch_id);
      unsorted_faces.erase(explored_facet);
      
      // Add/Remove facet's halfedge to/from explore list
      BOOST_FOREACH(halfedge_descriptor hd,
                    halfedges_around_face(halfedge(explored_facet,polyhedron),
                                          polyhedron))
      {
        halfedge_descriptor current_he = hd;
        
        // do not explore heh again
        if ( current_he == he ) { continue; }
        
        // if current_he is not in to_explore set, add it, otherwise remove it
        // (because we just explore the facet he_begin is pointing to)
        if ( he_to_explore.erase(current_he) == 0 )
        {
          he_to_explore.insert(opposite(current_he,polyhedron));
        }
      }
    }
  }
}
 
} // end namespace PMP
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYHEDRA_H

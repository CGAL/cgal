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

#ifndef CGAL_MESH_3_DETECT_FEATURES_IN_POLYHEDRA_H
#define CGAL_MESH_3_DETECT_FEATURES_IN_POLYHEDRA_H

#include <CGAL/Mesh_3/Detect_features_in_polyhedra_fwd.h>
#include <set>

namespace CGAL {
namespace Mesh_3 {

template <typename Polyhedron>
void detect_features(Polyhedron& p,
                     typename Polyhedron::Traits::FT angle_in_deg)
{
  Detect_features_in_polyhedra<Polyhedron> go;
  go.detect_sharp_edges(p,angle_in_deg);
  go.detect_surface_patches(p);
  go.detect_vertices_incident_patches(p);
}

  
template <typename Polyhedron_>
class Detect_features_in_polyhedra
{
  typedef Polyhedron_ Polyhedron;
public:
  typedef typename Polyhedron::Traits       Geom_traits;
  typedef typename Geom_traits::Vector_3    Vector_3;
  typedef typename Geom_traits::FT          FT;
  
  typedef typename Polyhedron::Halfedge_handle  Halfedge_handle;
  typedef typename Polyhedron::Facet_handle     Facet_handle;
  typedef typename Polyhedron::Halfedge         Halfedge;
  typedef typename Polyhedron::Facet            Facet;
  
  typedef std::set<Facet*>      Facet_handle_set;
  typedef std::set<Halfedge*>   He_handle_set;
  
public:
  Detect_features_in_polyhedra() : current_surface_index_(1) {}
  
  void detect_sharp_edges(Polyhedron& polyhedron,
                          FT angle_in_deg = FT(60)) const;
  void detect_surface_patches(Polyhedron& polyhedron);
  void detect_vertices_incident_patches(Polyhedron& p);
  
private:
  Vector_3 facet_normal(const Facet_handle& f) const;
  bool is_sharp(const Halfedge_handle& he, FT cos_angle) const;
  void flood(Facet& f, const int index,
             Facet_handle_set& unsorted_faces) const;
  
private:
  // Stores the current surface index (usefull to detect different patches
  // on different polyhedra)
  int current_surface_index_;
};

  
template <typename P_>
void
Detect_features_in_polyhedra<P_>::
detect_sharp_edges(Polyhedron& polyhedron, FT angle_in_deg) const
{
  // Initialize vertices
  for(typename Polyhedron::Vertex_iterator v = polyhedron.vertices_begin(),
      end = polyhedron.vertices_end() ; v != end ; ++v)
  {
    v->nb_of_feature_edges = 0;
  }
  
  FT cos_angle ( std::cos(CGAL::to_double(angle_in_deg) * CGAL_PI / 180.) );
  
  // Detect sharp edges
  for(typename Polyhedron::Halfedge_iterator he = polyhedron.edges_begin(),
      end = polyhedron.edges_end() ; he != end ; ++he)
  {
    if(he->is_border() || is_sharp(he,cos_angle))
    {
      he->set_feature_edge(true);
      he->opposite()->set_feature_edge(true);
      
      ++he->vertex()->nb_of_feature_edges;
      ++he->opposite()->vertex()->nb_of_feature_edges;
    }
  }
}


template <typename P_>
void
Detect_features_in_polyhedra<P_>::
detect_surface_patches(Polyhedron& polyhedron)
{
  // Initialize unsorted_faces
  Facet_handle_set unsorted_faces;
  for ( typename Polyhedron::Facet_iterator fit = polyhedron.facets_begin(),
       end = polyhedron.facets_end() ; fit != end ; ++fit )
  {
    unsorted_faces.insert(&*fit);
  }
  
  // Flood
  while ( ! unsorted_faces.empty() )
  {
    Facet& f = **(unsorted_faces.begin());
    unsorted_faces.erase(unsorted_faces.begin());
    
    f.set_patch_id(current_surface_index_);
    flood(f,current_surface_index_,unsorted_faces);
    ++current_surface_index_;
  }
}
  
  
template <typename P_>
void
Detect_features_in_polyhedra<P_>::
detect_vertices_incident_patches(Polyhedron& polyhedron)
{
  for( typename Polyhedron::Vertex_iterator vit = polyhedron.vertices_begin(),
      vend = polyhedron.vertices_end() ; vit != vend ; ++vit )
  {
    // Look only at feature vertices
    if( ! vit->is_feature_vertex() ) { continue; }
    
    // Loop on incident facets of vit
    typename Polyhedron::Halfedge_around_vertex_const_circulator
      he = vit->vertex_begin(), he_end(he);
    do
    {
      if( ! he->is_border() )
      {
        vit->add_incident_patch(he->facet()->patch_id());
      }
      else if( ! he->opposite()->is_border() )
      {
        vit->add_incident_patch(he->opposite()->facet()->patch_id());
      }
    } while(++he != he_end);
  }
}
  
// -----------------------------------
// Private methods
// -----------------------------------
template <typename P_>
typename Detect_features_in_polyhedra<P_>::Vector_3
Detect_features_in_polyhedra<P_>::
facet_normal(const Facet_handle& f) const
{
  Vector_3 sum = CGAL::NULL_VECTOR;
  typename Facet::Halfedge_around_facet_circulator h = f->facet_begin();
  
  do
  {
    Vector_3 normal = CGAL::cross_product(
      h->next()->vertex()->point() - h->vertex()->point(), 
      h->next()->next()->vertex()->point() - h->next()->vertex()->point());
    
    FT sqnorm = normal * normal;
    if ( ! CGAL_NTS is_zero(sqnorm) )
    {
      normal = normal / CGAL::sqrt(sqnorm);
      sum = sum + normal;
    }
  }
  while (++h != f->facet_begin());
  
  FT sqnorm = sum * sum;
  
  return (! CGAL_NTS is_zero(sqnorm)) ? sum / CGAL::sqrt(sqnorm)
                                      : CGAL::NULL_VECTOR;
}


template <typename P_>
bool
Detect_features_in_polyhedra<P_>::
is_sharp(const Halfedge_handle& he, FT cos_angle) const
{
  Facet_handle f1 = he->facet();
  Facet_handle f2 = he->opposite()->facet();
  if(f1 == NULL || f2 == NULL)
    return false;
  
  const Vector_3& n1 = facet_normal(f1);
  const Vector_3& n2 = facet_normal(f2);
  
  if ( n1 * n2 <= cos_angle )
    return true;
  else
    return false;
}
  
 
template <typename P_>
void
Detect_features_in_polyhedra<P_>::
flood(Facet& f, const int index, Facet_handle_set& unsorted_faces) const
{
  typedef typename Facet::Halfedge_around_facet_circulator Facet_he_circ;
  
  Facet_he_circ begin = f.facet_begin();
  Facet_he_circ done = begin;
  
  // Initialize he_to_explore with halfedges of the starting facet
  He_handle_set he_to_explore;
  CGAL_For_all(begin,done)
  {
    he_to_explore.insert(&*(begin->opposite()));
  }
  
  // While there is something to explore
  while ( ! he_to_explore.empty() )
  {
    // Get next halfedge to explore
    Halfedge& he = **(he_to_explore.begin());
    he_to_explore.erase(he_to_explore.begin());
    
    // If we don't go through a border of the patch
    if ( ! he.is_feature_edge() && ! he.is_border() )
    {
      Facet& explored_facet = *(he.facet());
      
      // Mark facet and delete it from unsorted
      explored_facet.set_patch_id(index);
      unsorted_faces.erase(&explored_facet);
      
      // Add/Remove facet's halfedge to/from explore list
      Facet_he_circ he_begin = explored_facet.facet_begin();
      Facet_he_circ he_done = he_begin;
      
      CGAL_For_all(he_begin,he_done)
      {
        Halfedge& current_he = *he_begin;
        
        // do not explore heh again
        if ( &current_he == &he ) { continue; }
        
        // if current_he is not in to_explore set, add it, otherwise remove it
        // (because we just explore the facet he_begin is pointing to)
        if ( he_to_explore.erase(&current_he) == 0 )
        {
          he_to_explore.insert(&*(current_he.opposite()));
        }
      }
    }
  }
}
 
} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_DETECT_FEATURES_IN_POLYHEDRA_H

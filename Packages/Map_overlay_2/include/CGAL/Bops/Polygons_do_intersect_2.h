// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Eti Ezra          <estere@post.tau.ac.il>
#ifndef CGAL_POLYGONS_DO_INTERSECT_2_H
#define CGAL_POLYGONS_DO_INTERSECT_2_H

CGAL_BEGIN_NAMESPACE

template <class Traits_>
class Polygons_do_intersect_2
{
  typedef Traits_       Traits;
  typedef typename Traits::Point_2           Point_2;
  typedef typename Traits::X_monotone_curve_2         X_monotone_curve_2;
  typedef typename Traits::Curve_2           Curve_2;
  typedef Holes_split_dcel<Traits>           Dcel;
  typedef Planar_map_2<Dcel,Traits>          Planar_map;
  typedef Map_overlay_2<Planar_map>          MapOverlay;
  typedef Boolean_operations_2<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
public:
  
  template <class Polygon> 
  bool operator()(const Polygon& polygon1,
                  const Polygon& polygon2) const
  { 
    PmWalkPL pm_walk1, pm_walk2;
    Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
    
    std::list<Curve_2> edges1, edges2;
    
    typename Polygon::Edge_const_iterator e_iter;
    for (e_iter = polygon1.edges_begin(); 
         e_iter != polygon1.edges_end(); ++e_iter)
      edges1.push_back(*e_iter);
    //edges1.push_back(traits.tsegment(*e_iter));
    
    for (e_iter = polygon2.edges_begin(); 
         e_iter != polygon2.edges_end(); ++e_iter)
      edges2.push_back(*e_iter);
    //edges2.push_back(traits.segment(*e_iter));
    
    pm1.insert(edges1.begin(), edges1.end());
    pm2.insert(edges2.begin(), edges2.end());
    
    pm1.unbounded_face()->set_ignore_bop(false); 
    pm2.unbounded_face()->set_ignore_bop(false);

    PmWalkPL   ovl_walk;
    //   MapOverlay map1(pm1);
    // MapOverlay map2(pm2);
    Bops bops(MapOverlay(pm1, pm2, &ovl_walk));

    Faces_container     faces;
    Halfedges_container halfedges;
    Vertices_container  vertices;
    
    bops.intersection(faces,halfedges,vertices);
    
    return (!faces.empty() ||
            !halfedges.empty() ||
            !vertices.empty()) ;
  }
};

CGAL_END_NAMESPACE

#endif







// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-11 $
// release_date  : $CGAL_Date: 2002/08/04 $
//
// file          : include/CGAL/Polygons_do_intersect_2.h
// package       : Map_overlay (1.12)
// maintainer    : Efi Fogel <efif@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra          <estere@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_POLYGON_DO_INTERSECT_2_H
#define CGAL_POLYGON_DO_INTERSECT_2_H

#include <CGAL/Pm_walk_along_line_point_location.h>

#ifndef PLANAR_MAP_2
#include <CGAL/Planar_map_2.h>
#endif

#ifndef CGAL_BOP_DEFAULT_DCEL_H
#include <CGAL/Bop_default_dcel.h>
#endif

#ifndef CGAL_MAP_OVERLAY_DEFAULT_NOTIFIER_H
#include <CGAL/Map_overlay_default_notifier.h>
#endif

#ifndef CGAL_MAP_OVERLAY_H
#include <CGAL/Map_overlay.h>
#endif

#ifndef BOOLEAN_OPERATIONS_2_H
#include <CGAL/Boolean_operations_2.h>
#endif

#include <CGAL/sweep_to_construct_planar_map_2.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Bops/Holes_split_dcel.h>
#include <CGAL/Bops/Holes_split_notifier.h>
#include <CGAL/Bops/Holes_split.h>
#include <CGAL/Bops/Polygons_from_faces.h>

CGAL_BEGIN_NAMESPACE

template <class Traits_>
class Polygons_do_intersect_2
{
  typedef Traits_       Traits;
  typedef typename Traits::Point_2           Point_2;
  typedef typename Traits::X_curve_2         X_curve_2;
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
  void operator()(const Polygon& polygon1,
                  const Polygon& polygon2) const
  { 
    PmWalkPL pm_walk1, pm_walk2;
    Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
    
    Traits traits;
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




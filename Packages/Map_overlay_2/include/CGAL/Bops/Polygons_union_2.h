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
// file          : include/CGAL/Polygons_union_2.h
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
#ifndef CGAL_POLYGONS_UNION_2_H
#define CGAL_POLYGONS_UNION_2_H

CGAL_BEGIN_NAMESPACE

template <class Traits_>
class Polygons_union_2
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
  
  template <class Polygon,
    class PolygonOutputIterator, 
    class CurvesOutputIterator,
    class PointsOutputIterator> 
  void operator()(const Polygon& polygon1,
                  const Polygon& polygon2,
                  PolygonOutputIterator polygons, 
                  CurvesOutputIterator curves,
                  PointsOutputIterator points) const
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
    
    bops.Union(faces,halfedges,vertices);
  
    typename Vertices_container::iterator v_iter;
    for(v_iter = vertices.begin(); 
        v_iter != vertices.end(); ++v_iter, ++points)
      *points = (*v_iter)->point();

    typename Halfedges_container::iterator h_iter;
    for(h_iter = halfedges.begin(); 
        h_iter != halfedges.end(); ++h_iter, ++h_iter, ++curves)
      *curves = (*h_iter)->curve();
    
    std::list<typename Traits::Curve> edges;
    typename Faces_container::iterator f_iter;
    for(f_iter = faces.begin(); f_iter != faces.end(); ++f_iter)
      {
        if ((*f_iter)->is_unbounded())
          continue;
      
        typename Planar_map::Holes_const_iterator hit;
        for (hit = (*f_iter)->holes_begin(); 
             hit != (*f_iter)->holes_end(); ++hit) 
          {
            typename  Planar_map::Ccb_halfedge_const_circulator cc(*hit);
            do {
              edges.push_back(cc->curve());
              //cout<<"hole:"<<cc->curve()<<endl;
            } while (++cc != *hit);
          }
        
        typename Planar_map::Ccb_halfedge_const_circulator cc = 
          (*f_iter)->outer_ccb();
        do {
          edges.push_back(cc->curve());
          //cout<<"outer ccb: ("<< CGAL::to_double(cc->curve().source().x()) <<","
          //    << CGAL::to_double(cc->curve().source().y()) <<") ("
          //    << CGAL::to_double(cc->curve().target().x()) << ","
          //    << CGAL::to_double(cc->curve().target().y()) << ")" <<endl;
        } while (++cc != (*f_iter)->outer_ccb());
        
        // Generating Polygons:
        PmWalkPL pm_walk;
        Planar_map pm(&pm_walk);
        CGAL::sweep_to_construct_planar_map_2(edges.begin(),edges.end(), 
                                              traits, pm);
    
        // First, emanating vertical walls from each hole in the map.
        Notifier  notf;
        Holes_split<Planar_map,Notifier> holes_split(&traits);
        holes_split.split_holes(pm, notf);
    
        // Next, turning each face to a polygon.
        Polygons_from_faces<Planar_map, Polygon>  polygon_from_faces;
        polygon_from_faces(pm.faces_begin(), pm.faces_end(), polygons);
      }
  }
};

CGAL_END_NAMESPACE

#endif







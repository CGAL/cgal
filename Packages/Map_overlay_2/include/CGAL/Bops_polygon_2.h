#ifndef CGAL_BOPS_POLYGON_2_H
#define CGAL_BOPS_POLYGON_2_H

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

template < class Polygon, class ArrangementTraits_2 >
bool do_intersect(const Polygon& A,
                  const Polygon& B, 
                  ArrangementTraits_2& traits)
{
  typedef ArrangementTraits_2               Traits;
  typedef Bop_default_dcel<Traits>          Dcel;
  typedef Planar_map_2<Dcel,Traits>         Planar_map;
  typedef Map_overlay_2<Planar_map>         MapOverlay;
  typedef Boolean_operations_2<MapOverlay>  Bops;
  typedef CGAL::Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  
  typedef Bops::Faces_container           Faces_container;
  typedef Bops::Halfedges_container       Halfedges_container;
  typedef Bops::Vertices_container        Vertices_container;

  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_itertor edge_iter;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter)
    pm1.insert(*edge_iter);
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter)
    pm2.insert(*edge_iter);
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bop.intersection(faces,halfedges,vertices);
  
  return (!faces.empty() ||
          !halfedges.empty() ||
          !vertices.empty()) ;
}

template < class Polygon, 
           class ArrangementTraits_2, 
           class BopsTraits_2,
           class PolygonOutputIterator, 
           class CurvesOutputIterator,
           class PointsOutputIterator >
void intersection(const Polygon& A,
                  const Polygon& B, 
                  ArrangementTraits_2& traits,
                  BopsTraits_2&  bops_traits,
                  PolygonOutputIterator polygons, 
                  CurvesOutputIterator curves,
                  PointsOutputIterator points)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations_2<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  pm1.insert(bops_traits.edges_begin(A), bops_traits.edges_end(A));
  pm2.insert(bops_traits.edges_begin(B), bops_traits.edges_end(B));

  /*typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
    }*/
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.intersection(faces,halfedges,vertices);
  
  for(typename Vertices_container::iterator v_iter = vertices.begin(); 
      v_iter != vertices.end(); ++v_iter, ++points)
    *points = (*v_iter)->point();


  for(typename Halfedges_container::iterator h_iter = halfedges.begin(); 
      h_iter != halfedges.end(); ++h_iter, ++h_iter, ++curves)
    *curves = (*h_iter)->curve();
  
  std::list<typename Traits::Curve> edges;
  for(typename Faces_container::iterator f_iter = faces.begin(); 
      f_iter != faces.end(); ++f_iter){
    
    if ((*f_iter)->is_unbounded())
      continue;
     
    for (typename  Planar_map::Holes_const_iterator hit = (*f_iter)->holes_begin(); 
         hit != (*f_iter)->holes_end(); ++hit) {
      typename  Planar_map::Ccb_halfedge_const_circulator cc(*hit);
      do {
        edges.push_back(cc->curve());
        cout<<"hole:"<<cc->curve()<<endl;
      } while (++cc != *hit);
    }
    
    typename Planar_map::Ccb_halfedge_const_circulator cc=(*f_iter)->outer_ccb();
    do {
      edges.push_back(cc->curve());
      cout<<"outer ccb: ("<< CGAL::to_double(cc->curve().source().x()) <<","
          << CGAL::to_double(cc->curve().source().y()) <<") ("
          << CGAL::to_double(cc->curve().target().x()) << ","
          << CGAL::to_double(cc->curve().target().y()) << ")" <<endl;
    } while (++cc != (*f_iter)->outer_ccb());
    
    // Generating Polygons:
    PmWalkPL pm_walk;
    Planar_map pm(&pm_walk);
    CGAL::sweep_to_construct_planar_map_2(edges.begin(),edges.end(), 
                                          traits, pm);
    
    // First, aminating vertical walls from each hole in the map.
    Notifier  notf;
    Holes_split<Planar_map,Notifier,BopsTraits_2> holes_split(&bops_traits);
    holes_split.split_holes(pm, notf);
    
    // Next, turning each face to a polygon.
    Polygons_from_faces<Planar_map, Polygon>  polygon_from_faces;
    polygon_from_faces(pm.faces_begin(), pm.faces_end(), polygons);
  }
}

/*
  template < class Polygon, class ArrangementTraits_2, 
           class PolygonOutputIterator, 
           class CurvesOutputIterator>
void intersection(const Polygon& A,
                  const Polygon& B, 
                  ArrangementTraits_2& traits,
                  PolygonOutputIterator polygons, 
                  CurvesOutputIterator curves)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
  }
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.intersection(faces,halfedges,vertices);

  for(typename Halfedges_container::iterator h_iter = halfedges.begin(); 
      h_iter != halfedges.end(); ++h_iter, ++h_iter, ++curves)
    *curves = (*h_iter)->curve();
  
  std::list<typename Traits::Curve> edges;
  for(typename Faces_container::iterator f_iter = faces.begin(); 
      f_iter != faces.end(); ++f_iter){
    
    if ((*f_iter)->is_unbounded())
      continue;
     
    for (typename  Planar_map::Holes_const_iterator hit = (*f_iter)->holes_begin(); 
         hit != (*f_iter)->holes_end(); ++hit) {
      typename  Planar_map::Ccb_halfedge_const_circulator cc(*hit);
      do {
        edges.push_back(cc->curve());
      } while (++cc != *hit);
    }
    
    typename Planar_map::Ccb_halfedge_const_circulator cc=(*f_iter)->outer_ccb();
    do {
      edges.push_back(cc->curve());
    } while (++cc != (*f_iter)->outer_ccb());
    
    // Generating Polygons:
    PmWalkPL pm_walk;
    Planar_map pm(&pm_walk);
    CGAL::sweep_to_construct_planar_map_2(edges.begin(),edges.end(), 
                                          traits, pm);
    
    // First, aminating vertical walls from each hole in the map.
    Notifier  notf;
    Holes_split<Planar_map,Notifier> holes_split;
    holes_split.split_holes(pm, notf);
    
    // Next, turning each face to a polygon.
    Polygons_from_faces<Planar_map, Polygon>  polygon_from_faces;
    polygon_from_faces(pm.faces_begin(), pm.faces_end(), polygons);
  }
}

template < class Polygon, class ArrangementTraits_2, 
           class PolygonOutputIterator, 
           class CurvesOutputIterator,
           class PointsOutputIterator >
void intersection(const Polygon& A,
                  const Polygon& B, 
                  ArrangementTraits_2& traits,
                  PolygonOutputIterator polygons, 
                  PointsOutputIterator points)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
  }
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.intersection(faces,halfedges,vertices);
  
  for(typename Vertices_container::iterator v_iter = vertices.begin(); 
      v_iter != vertices.end(); ++v_iter, ++points)
    *points = (*v_iter)->point();
  
  std::list<typename Traits::Curve> edges;
  for(typename Faces_container::iterator f_iter = faces.begin(); 
      f_iter != faces.end(); ++f_iter){
    
    if ((*f_iter)->is_unbounded())
      continue;
     
    for (typename  Planar_map::Holes_const_iterator hit = (*f_iter)->holes_begin(); 
         hit != (*f_iter)->holes_end(); ++hit) {
      typename  Planar_map::Ccb_halfedge_const_circulator cc(*hit);
      do {
        edges.push_back(cc->curve());
      } while (++cc != *hit);
    }
    
    typename Planar_map::Ccb_halfedge_const_circulator cc=(*f_iter)->outer_ccb();
    do {
      edges.push_back(cc->curve());
    } while (++cc != (*f_iter)->outer_ccb());
    
    // Generating Polygons:
    PmWalkPL pm_walk;
    Planar_map pm(&pm_walk);
    CGAL::sweep_to_construct_planar_map_2(edges.begin(),edges.end(), 
                                          traits, pm);
    
    // First, aminating vertical walls from each hole in the map.
    Notifier  notf;
    Holes_split<Planar_map,Notifier> holes_split;
    holes_split.split_holes(pm, notf);
    
    // Next, turning each face to a polygon.
    Polygons_from_faces<Planar_map, Polygon>  polygon_from_faces;
    polygon_from_faces(pm.faces_begin(), pm.faces_end(), polygons);
  }
}

template < class Polygon, class ArrangementTraits_2, 
           class PolygonOutputIterator>
void intersection(const Polygon& A,
                  const Polygon& B, 
                  ArrangementTraits_2& traits,
                  PolygonOutputIterator polygons)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
  }
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.intersection(faces,halfedges,vertices);
  
  std::list<typename Traits::Curve> edges;
  for(typename Faces_container::iterator f_iter = faces.begin(); 
      f_iter != faces.end(); ++f_iter){
    
    if ((*f_iter)->is_unbounded())
      continue;
     
    for (typename  Planar_map::Holes_const_iterator hit = (*f_iter)->holes_begin(); 
         hit != (*f_iter)->holes_end(); ++hit) {
      typename  Planar_map::Ccb_halfedge_const_circulator cc(*hit);
      do {
        edges.push_back(cc->curve());
      } while (++cc != *hit);
    }
    
    typename Planar_map::Ccb_halfedge_const_circulator cc=(*f_iter)->outer_ccb();
    do {
      edges.push_back(cc->curve());
    } while (++cc != (*f_iter)->outer_ccb());
    
    // Generating Polygons:
    PmWalkPL pm_walk;
    Planar_map pm(&pm_walk);
    CGAL::sweep_to_construct_planar_map_2(edges.begin(),edges.end(), 
                                          traits, pm);
    
    // First, aminating vertical walls from each hole in the map.
    Notifier  notf;
    Holes_split<Planar_map,Notifier> holes_split;
    holes_split.split_holes(pm, notf);
    
    // Next, turning each face to a polygon.
    Polygons_from_faces<Planar_map, Polygon>  polygon_from_faces;
    polygon_from_faces(pm.faces_begin(), pm.faces_end(), polygons);
  }
}

template < class Polygon, class ArrangementTraits_2, 
           class CurvesOutputIterator,
           class PointsOutputIterator >
void intersection(const Polygon& A,
                  const Polygon& B, 
                  ArrangementTraits_2& traits,
                  CurvesOutputIterator curves,
                  PointsOutputIterator points)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
  }
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.intersection(faces,halfedges,vertices);
  
  for(typename Vertices_container::iterator v_iter = vertices.begin(); 
      v_iter != vertices.end(); ++v_iter, ++points)
    *points = (*v_iter)->point();


  for(typename Halfedges_container::iterator h_iter = halfedges.begin(); 
      h_iter != halfedges.end(); ++h_iter, ++h_iter, ++curves)
    *curves = (*h_iter)->curve();
}

template < class Polygon, class ArrangementTraits_2, 
           class CurvesOutputIterator>
void intersection(const Polygon& A,
                  const Polygon& B, 
                  ArrangementTraits_2& traits,
                  CurvesOutputIterator curves)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
  }
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.intersection(faces,halfedges,vertices);

  for(typename Halfedges_container::iterator h_iter = halfedges.begin(); 
      h_iter != halfedges.end(); ++h_iter, ++h_iter, ++curves)
    *curves = (*h_iter)->curve();
}

template < class Polygon, class ArrangementTraits_2, 
           class PointsOutputIterator >
void intersection(const Polygon& A,
                  const Polygon& B, 
                  ArrangementTraits_2& traits,
                  PointsOutputIterator points)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
  }
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.intersection(faces,halfedges,vertices);
  
  for(typename Vertices_container::iterator v_iter = vertices.begin(); 
      v_iter != vertices.end(); ++v_iter, ++points)
    *points = (*v_iter)->point();
}*/

template < class Polygon, 
           class ArrangementTraits_2, 
           class BopsTraits_2,
           class PolygonOutputIterator, 
           class CurvesOutputIterator,
           class PointsOutputIterator >
void Union(const Polygon& A,
           const Polygon& B, 
           ArrangementTraits_2& traits,
           BopsTraits_2&  bops_traits,
           PolygonOutputIterator polygons, 
           CurvesOutputIterator curves,
           PointsOutputIterator points)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations_2<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
  }
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.Union(faces,halfedges,vertices);
  
  for(typename Vertices_container::iterator v_iter = vertices.begin(); 
      v_iter != vertices.end(); ++v_iter, ++points)
    *points = (*v_iter)->point();


  for(typename Halfedges_container::iterator h_iter = halfedges.begin(); 
      h_iter != halfedges.end(); ++h_iter, ++h_iter, ++curves)
    *curves = (*h_iter)->curve();
  
  std::list<typename Traits::Curve> edges;
  for(typename Faces_container::iterator f_iter = faces.begin(); 
      f_iter != faces.end(); ++f_iter){
    
    if ((*f_iter)->is_unbounded())
      continue;
     
    for (typename  Planar_map::Holes_const_iterator hit = (*f_iter)->holes_begin(); 
         hit != (*f_iter)->holes_end(); ++hit) {
      typename  Planar_map::Ccb_halfedge_const_circulator cc(*hit);
      do {
        edges.push_back(cc->curve());
        cout<<"hole:"<<cc->curve()<<endl;
      } while (++cc != *hit);
    }
    
    typename Planar_map::Ccb_halfedge_const_circulator cc=(*f_iter)->outer_ccb();
    do {
      edges.push_back(cc->curve());
      cout<<"outer ccb: ("<< CGAL::to_double(cc->curve().source().x()) <<","
          << CGAL::to_double(cc->curve().source().y()) <<") ("
          << CGAL::to_double(cc->curve().target().x()) << ","
          << CGAL::to_double(cc->curve().target().y()) << ")" <<endl;
    } while (++cc != (*f_iter)->outer_ccb());
    
    // Generating Polygons:
    PmWalkPL pm_walk;
    Planar_map pm(&pm_walk);
    CGAL::sweep_to_construct_planar_map_2(edges.begin(),edges.end(), 
                                          traits, pm);
    
    // First, aminating vertical walls from each hole in the map.
    Notifier  notf;
    Holes_split<Planar_map,Notifier,BopsTraits_2> holes_split(&bops_traits);
    holes_split.split_holes(pm, notf);
    
    // Next, turning each face to a polygon.
    Polygons_from_faces<Planar_map, Polygon>  polygon_from_faces;
    polygon_from_faces(pm.faces_begin(), pm.faces_end(), polygons);
  }
}

template < class Polygon, 
           class ArrangementTraits_2, 
           class BopsTraits_2,
           class PolygonOutputIterator, 
           class CurvesOutputIterator,
           class PointsOutputIterator >
void symmetric_difference(const Polygon& A,
                          const Polygon& B, 
                          ArrangementTraits_2& traits,
                          BopsTraits_2&  bops_traits,
                          PolygonOutputIterator polygons, 
                          CurvesOutputIterator curves,
                          PointsOutputIterator points)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations_2<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
  }
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.symmetric_difference(faces,halfedges,vertices);
  
  for(typename Vertices_container::iterator v_iter = vertices.begin(); 
      v_iter != vertices.end(); ++v_iter, ++points)
    *points = (*v_iter)->point();


  for(typename Halfedges_container::iterator h_iter = halfedges.begin(); 
      h_iter != halfedges.end(); ++h_iter, ++h_iter, ++curves)
    *curves = (*h_iter)->curve();
  
  std::list<typename Traits::Curve> edges;
  for(typename Faces_container::iterator f_iter = faces.begin(); 
      f_iter != faces.end(); ++f_iter){
    
    if ((*f_iter)->is_unbounded())
      continue;
     
    for (typename  Planar_map::Holes_const_iterator hit = (*f_iter)->holes_begin(); 
         hit != (*f_iter)->holes_end(); ++hit) {
      typename  Planar_map::Ccb_halfedge_const_circulator cc(*hit);
      do {
        edges.push_back(cc->curve());
        cout<<"hole:"<<cc->curve()<<endl;
      } while (++cc != *hit);
    }
    
    typename Planar_map::Ccb_halfedge_const_circulator cc=(*f_iter)->outer_ccb();
    do {
      edges.push_back(cc->curve());
      cout<<"outer ccb: ("<< CGAL::to_double(cc->curve().source().x()) <<","
          << CGAL::to_double(cc->curve().source().y()) <<") ("
          << CGAL::to_double(cc->curve().target().x()) << ","
          << CGAL::to_double(cc->curve().target().y()) << ")" <<endl;
    } while (++cc != (*f_iter)->outer_ccb());
    
    // Generating Polygons:
    PmWalkPL pm_walk;
    Planar_map pm(&pm_walk);
    CGAL::sweep_to_construct_planar_map_2(edges.begin(),edges.end(), 
                                          traits, pm);
    
    // First, aminating vertical walls from each hole in the map.
    Notifier  notf;
    Holes_split<Planar_map,Notifier,BopsTraits_2> holes_split(&bops_traits);
    holes_split.split_holes(pm, notf);
    
    // Next, turning each face to a polygon.
    Polygons_from_faces<Planar_map, Polygon>  polygon_from_faces;
    polygon_from_faces(pm.faces_begin(), pm.faces_end(), polygons);
  }
}

template < class Polygon, 
           class ArrangementTraits_2, 
           class BopsTraits_2,
           class PolygonOutputIterator, 
           class CurvesOutputIterator,
           class PointsOutputIterator >
void difference(const Polygon& A,
                const Polygon& B, 
                ArrangementTraits_2& traits,
                BopsTraits_2&  bops_traits,
                PolygonOutputIterator polygons, 
                CurvesOutputIterator curves,
                PointsOutputIterator points)
{
  typedef ArrangementTraits_2              Traits;
  typedef Holes_split_dcel<Traits>         Dcel;
  typedef Planar_map_2<Dcel,Traits>        Planar_map;
  typedef Map_overlay<Planar_map>          MapOverlay;
  typedef Boolean_operations_2<MapOverlay>   Bops;
  typedef Pm_walk_along_line_point_location<Planar_map>   PmWalkPL;
  typedef Holes_split_notifier<Planar_map>         Notifier;
  
  typedef typename Bops::Faces_container           Faces_container;
  typedef typename Bops::Halfedges_container       Halfedges_container;
  typedef typename Bops::Vertices_container        Vertices_container;
  
  //Planar_map pm1(A.edges_begin(), A.edges_end());
  //Planar_map pm2(B.edges_begin(), B.edges_end());
  PmWalkPL pm_walk1, pm_walk2;
  Planar_map pm1(&pm_walk1), pm2(&pm_walk2);
  
  typename Polygon::Edge_const_iterator edge_iter;
  //  cout<<"first polygon"<<endl;
  for (edge_iter = A.edges_begin(); edge_iter != A.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm1.insert(*edge_iter);
  }
  //cout<<"second polygon"<<endl;
  for (edge_iter = B.edges_begin(); edge_iter != B.edges_end(); ++edge_iter){
    //  cout<<*edge_iter<<endl; //debug
    pm2.insert(*edge_iter);
  }
  
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);

  PmWalkPL   ovl_walk;
  MapOverlay map1(pm1);
  MapOverlay map2(pm2);
  Bops bops(MapOverlay(map1, map2, &ovl_walk));

  Faces_container     faces;
  Halfedges_container halfedges;
  Vertices_container  vertices;
  
  bops.difference(faces,halfedges,vertices);
  
  for(typename Vertices_container::iterator v_iter = vertices.begin(); 
      v_iter != vertices.end(); ++v_iter, ++points)
    *points = (*v_iter)->point();


  for(typename Halfedges_container::iterator h_iter = halfedges.begin(); 
      h_iter != halfedges.end(); ++h_iter, ++h_iter, ++curves)
    *curves = (*h_iter)->curve();
  
  std::list<typename Traits::Curve> edges;
  for(typename Faces_container::iterator f_iter = faces.begin(); 
      f_iter != faces.end(); ++f_iter){
    
    if ((*f_iter)->is_unbounded())
      continue;
     
    for (typename  Planar_map::Holes_const_iterator hit = (*f_iter)->holes_begin(); 
         hit != (*f_iter)->holes_end(); ++hit) {
      typename  Planar_map::Ccb_halfedge_const_circulator cc(*hit);
      do {
        edges.push_back(cc->curve());
        cout<<"hole:"<<cc->curve()<<endl;
      } while (++cc != *hit);
    }
    
    typename Planar_map::Ccb_halfedge_const_circulator cc=(*f_iter)->outer_ccb();
    do {
      edges.push_back(cc->curve());
      cout<<"outer ccb: ("<< CGAL::to_double(cc->curve().source().x()) <<","
          << CGAL::to_double(cc->curve().source().y()) <<") ("
          << CGAL::to_double(cc->curve().target().x()) << ","
          << CGAL::to_double(cc->curve().target().y()) << ")" <<endl;
    } while (++cc != (*f_iter)->outer_ccb());
    
    // Generating Polygons:
    PmWalkPL pm_walk;
    Planar_map pm(&pm_walk);
    CGAL::sweep_to_construct_planar_map_2(edges.begin(),edges.end(), 
                                          traits, pm);
    
    // First, aminating vertical walls from each hole in the map.
    Notifier  notf;
    Holes_split<Planar_map,Notifier,BopsTraits_2> holes_split(&bops_traits);
    holes_split.split_holes(pm, notf);
    
    // Next, turning each face to a polygon.
    Polygons_from_faces<Planar_map, Polygon>  polygon_from_faces;
    polygon_from_faces(pm.faces_begin(), pm.faces_end(), polygons);
  }
}

CGAL_END_NAMESPACE

//#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
//#include <CGAL/bops_Polygon_2.C>
//#endif

#endif // CGAL_BOPS_POLYGON_2_H

#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 400)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.0 or higher installed!\n";
 std::cout << "A LEDA version >= 4.0 is required to run GeoWin!\n";
 return 0;
}
#else 

#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h> 
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/leda_integer.h>

typedef leda_integer coord_type;
typedef CGAL::Cartesian<coord_type>  My_REP;

typedef CGAL::Point_2< My_REP >     My_CGALPoint;
typedef std::list<My_CGALPoint>     My_CGALPointlist;
typedef CGAL::Segment_2< My_REP >   My_CGALSegment;
typedef std::list<My_CGALSegment>   My_CGALSegmentlist;
typedef CGAL::Ray_2< My_REP >       My_CGALRay;
typedef std::list<My_CGALRay>       My_CGALRaylist;

typedef CGAL::Triangulation_euclidean_traits_2<My_REP> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>  Triangulation_2;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Delaunay_triangulation_2;

typedef Triangulation_2::Edge Edge;
typedef Triangulation_2::Locate_type Locate_type;
typedef Triangulation_2::Edge_iterator  Edge_iterator;

#include <CGAL/Object.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/geowin_support.h>

#if defined GEOWIN_SUPPORT_NO_TEMPLATES
#include "geowin_workaround1.h"
#endif

class geo_triang : public geowin_update<std::list<My_CGALPoint>, std::list<My_CGALSegment> >
{
public:
 void update(const My_CGALPointlist& L, My_CGALSegmentlist& Sl)
 {
  Triangulation_2 tr;    
  Sl.clear();      
                   
  std::list<My_CGALPoint>::const_iterator it;
  it= L.begin();
  My_CGALPoint pakt;
 
  for (; it != L.end() ; ++it) {
        pakt= *it;
        tr.push_back(pakt);                            
  }

  Edge_iterator eit = tr.edges_begin();
  Edge_iterator beyond = tr.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;         
       Sl.push_back(tr.segment(eact));               
       ++eit;  
  }         
 }
};

class geo_delaunay_triang : public geowin_update<std::list<My_CGALPoint>, std::list<My_CGALSegment> >
{
public:
 void update(const My_CGALPointlist& L, My_CGALSegmentlist& Sl)
 {
  Delaunay_triangulation_2 dt;    
  Sl.clear();      

  dt.insert(L.begin(),L.end());

  Edge_iterator eit = dt.edges_begin();
  Edge_iterator beyond = dt.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;         
       Sl.push_back(dt.segment(eact));               
       ++eit;  
  }       
 }
};

class geo_voro1 : public geowin_update<std::list<My_CGALPoint>, std::list<My_CGALSegment> >
{
public:
 void update(const My_CGALPointlist& L, My_CGALSegmentlist& Sl)
 {
  Delaunay_triangulation_2 dt;    
  Sl.clear();      
  dt.insert(L.begin(),L.end()); 

  Edge_iterator eit, est=dt.edges_begin(), eend=dt.edges_end();
  for (eit=est; eit != eend; ++eit){
    My_CGALSegment s;
    CGAL::Object o = dt.dual(eit);
    if (CGAL::assign(s,o)) Sl.push_back(s);
  }
 }
};

class geo_voro2 : public geowin_update<std::list<My_CGALPoint>, std::list<My_CGALRay> >
{
public:
 void update(const My_CGALPointlist& L, My_CGALRaylist& Sl)
 {
  Delaunay_triangulation_2 dt;    
  Sl.clear();      
  dt.insert(L.begin(),L.end()); 

  Edge_iterator eit, est=dt.edges_begin(), eend=dt.edges_end();
  for (eit=est; eit != eend; ++eit){
    My_CGALRay s;
    CGAL::Object o = dt.dual(eit);
    if (CGAL::assign(s,o)) Sl.push_back(s);
  }
 }
};

int main()
{
  geowin_init_default_type((My_CGALPointlist*)0, leda_string("My_CGALPointList"));

  My_CGALPointlist L;

  GeoWin GW("CGAL - Triangulation demo");
 
  // build new scene...
  geo_scene my_scene= GW.new_scene(L);  

  // build result scenes ...
  geo_triang triangulate;
  geo_scene result  = GW.new_scene(triangulate ,my_scene , leda_string("Triangulation"));
  GW.set_visible(result,true);
 
  geo_delaunay_triang deltria;
  geo_scene res2 = GW.new_scene(deltria, my_scene, leda_string("Delaunay Triangulation"));
  GW.set_color(res2, leda_red);
  GW.set_visible(res2, true);

  geo_voro1 voro_edges;
  geo_scene res3 = GW.new_scene(voro_edges, my_scene, leda_string("Voronoi Edges"));
  GW.set_color(res3, leda_blue);
  GW.set_visible(res3, true);

  geo_voro2 voro_rays;
  geo_scene res4 = GW.new_scene(voro_rays, my_scene, leda_string("Voronoi Edges"));
  GW.set_color(res4, leda_blue);
  GW.set_visible(res4, true);
  
  GW.edit(my_scene);
  
  return 0;  
}

#endif

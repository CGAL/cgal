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
#include <CGAL/geowin_support.h>

typedef double coord_type;
typedef CGAL::Cartesian<coord_type>  Rep;
typedef CGAL::Point_2<Rep>  Point;
typedef CGAL::Segment_2<Rep>  Segment;

typedef CGAL::Triangulation_euclidean_traits_2<Rep> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>  Triang_2;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Delaunay_triang_2;

typedef Delaunay_triang_2::Face::Face_handle Face_handle;
typedef Delaunay_triang_2::Line_face_circulator Line_face_circulator;

typedef Triang_2::Face  Face;
typedef Triang_2::Vertex Vertex;
typedef Triang_2::Edge Edge;
typedef Triang_2::Vertex_handle Vertex_handle;
typedef Triang_2::Edge_iterator  Edge_iterator;

class geo_delau : public geowin_update<std::list<CGALPoint>,std::list<CGALSegment> >
{
public:
 void update(const CGALPointlist& L, CGALSegmentlist& Sl)
 {
  Delaunay_triang_2 dt;    
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

class geo_triangles : public geowin_redraw, public geowin_update<CGALPointlist, CGALPointlist >
{
public:

  CGALTrianglelist LT;
  geo_scene lines;

  virtual ~geo_triangles() {}

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {  
    std::list<CGALTriangle>::const_iterator it   = LT.begin(), stop = LT.end();
    leda_color old = W.set_fill_color(leda_green);
    while( it != stop )
    {
       W << convert_to_leda(*it);
       it++;
    }
    W.set_fill_color(old);
  }
  
  virtual bool write_postscript(ps_file& PS,leda_color c1,leda_color c2)
  {
    std::list<CGALTriangle>::const_iterator it   = LT.begin(), stop = LT.end();
    leda_color old = PS.set_fill_color(leda_green);
    while( it != stop )
    {
       PS << convert_to_leda(*it); 
       it++;
    }
    PS.set_fill_color(old);  
    return false;
  }  

  virtual void update(const CGALPointlist& L, CGALPointlist&)
  { 

    Delaunay_triang_2 dt;    
    dt.insert(L.begin(),L.end());
    LT.clear();

    if (dt.dimension() != 2) return;

    GeoWin* gw = get_geowin(lines);

    CGALLinelist LST;
    gw->get_objects(lines,LST);

    std::list<CGALLine>::const_iterator it;
    CGALLine lakt;
    CGALPoint p1,p2; 

    for(it=LST.begin(); it != LST.end(); ++it) { 
       lakt= *it;
       p1=lakt.point(1); p2=lakt.point(2);
       Face_handle f = dt.locate(p1);
      
       Line_face_circulator lfc=dt.line_walk(p1,p2,f), done(lfc);
 
       if (lfc== (CGAL_NULL_TYPE) NULL) continue;

       do {
	  if(! dt.is_infinite( lfc  )){ LT.push_back(dt.triangle(lfc)); }
       }  
       while (++lfc != done); 
      
    }   
  }

};

class geo_triangles2 : public geowin_redraw, public geowin_update<CGALPointlist, CGALPointlist >
{
public:

  CGALTrianglelist LT;
  geo_scene locate_points;

  virtual ~geo_triangles2() {}

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {  
    std::list<CGALTriangle>::const_iterator it   = LT.begin(), stop = LT.end();
    leda_color old = W.set_fill_color(leda_red);
    while( it != stop )
    {
       W << convert_to_leda(*it);
       it++;
    }
    W.set_fill_color(old);
  }
  
  virtual bool write_postscript(ps_file& PS,leda_color c1,leda_color c2)
  {
    std::list<CGALTriangle>::const_iterator it   = LT.begin(), stop = LT.end();
    leda_color old = PS.set_fill_color(leda_red);
    while( it != stop )
    {
       PS << convert_to_leda(*it); 
       it++;
    }
    PS.set_fill_color(old);  
    return false;
  }  

  virtual void update(const CGALPointlist& L, CGALPointlist&)
  { 

    Delaunay_triang_2 dt;    
    dt.insert(L.begin(),L.end());

    LT.clear();

    if (dt.dimension() != 2) return;

    GeoWin* gw = get_geowin(locate_points);   

    CGALPointlist LST;
    gw->get_objects(locate_points,LST);

    std::list<CGALPoint>::const_iterator it;
    CGALPoint lakt;

    for(it=LST.begin(); it != LST.end(); ++it) { 
       lakt= *it;
       Face_handle f = dt.locate(lakt);
       if (f != NULL && !dt.is_infinite(f)) LT.push_back(dt.triangle(f));
    }   
  }
};


int main()
{
  geo_scene_traits<std::list<CGALPoint> > TR(geowin_info_fcn, geowin_IntersectsBox, geowin_BoundingBox, 
             geowin_Translate, geowin_Rotate, geowin_generate_objects, leda_string("CGALPointlist"));
  geo_scene_traits<std::list<CGALLine> > TR2(geowin_info_fcn, geowin_IntersectsBox, geowin_BoundingBox, 
             geowin_Translate, geowin_Rotate, geowin_generate_objects, leda_string("CGALLinelist"));	     
  geowin_init_default_type(TR);
  geowin_init_default_type(TR2);

  CGALPointlist L, LOC;
  CGALLinelist  CGLL;

  GeoWin GW("CGAL - Triangulation demo");
  geo_scene my_scene= GW.new_scene(L); 
  // another input scene for the lines ...
  geo_scene line_scene= GW.new_scene(CGLL);
  GW.set_color(line_scene, leda_grey2);

  // another input scene for the points for location algor. ...
  geo_scene pointloc_scene= GW.new_scene(LOC);
  GW.set_color(pointloc_scene, leda_blue);

  geo_triangles TRS;
  geo_scene sc1 = GW.new_scene( TRS, TRS, my_scene, "Triangles"); 
  GW.set_color(sc1, leda_blue);
  GW.set_visible(sc1,true);
  TRS.lines = line_scene;
 
  geo_delau delaunay_triang;
  geo_scene res2 = GW.new_scene(delaunay_triang, my_scene, leda_string("Delaunay Triangulation"));
  GW.set_color(res2, leda_blue);
  GW.set_line_width(res2, 2);

  geo_triangles2 TPT;
  geo_scene res3 = GW.new_scene( TPT, TPT, my_scene, "Triangles2"); 
  GW.set_color(res3, leda_red);
  TPT.locate_points = pointloc_scene;
  
  GW.add_dependence(line_scene,sc1);
  GW.add_dependence(pointloc_scene,res3);
  GW.set_all_visible(true);

  GW.edit(my_scene);
  return 0;  
}

#endif

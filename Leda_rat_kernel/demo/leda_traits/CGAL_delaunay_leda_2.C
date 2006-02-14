
#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 430)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.3 or higher installed!\n";
 std::cout << "A LEDA version >= 4.3 is required !\n";
 return 0;
}
#else 

#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/geowin_support.h>


#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits                           K;
typedef K::Point_2                                             Point;
typedef K::Segment_2                                           Segment;
typedef K::Line_2                                              Line;
typedef K::Triangle_2                                          Triangle;

typedef CGAL::Triangulation_2<K>                               Triang_2;
typedef CGAL::Delaunay_triangulation_2<K>                      Delaunay_triang_2;
typedef Delaunay_triang_2::Face::Face_handle                   Face_handle;
typedef Delaunay_triang_2::Line_face_circulator                Line_face_circulator;

typedef Triang_2::Edge                                         Edge;
typedef Triang_2::Vertex_handle                                Vertex_handle;
typedef Triang_2::Edge_iterator                                Edge_iterator;

Delaunay_triang_2 dt;


namespace CGAL {
leda_segment convert_to_leda(const leda_rat_segment& s) { return s.to_float(); }
leda_triangle convert_to_leda(const leda_rat_triangle& t) { return t.to_float(); }
}


class geo_delau : public geowin_update<std::list<Point>,std::list<Segment> > ,
                  public geowin_redraw
{
public:
 // support for incremental operations ...
 bool insert(const Point& p)
 {
  std::cout << "insert:" << p << "\n";
  dt.insert(p);
  return true; 
 }
 
 bool del(const Point& p)
 {
  std::cout << "del:" << p << "\n";
  Delaunay_triang_2::Vertex_handle vh = dt.nearest_vertex(p);
  dt.remove(vh);
  return true;  
 } 

 void draw(leda_window& W,leda_color c1,leda_color c2,double x1,double y1,double x2,double y2)
 {
  Edge_iterator eit = dt.edges_begin();
  Edge_iterator beyond = dt.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;   
       W.draw_segment(CGAL::convert_to_leda(dt.segment(eact)),c1);                    
       ++eit;  
  }    
 }

 void update(const std::list<Point>& L, std::list<Segment>&)
 {
  dt.clear();    
  dt.insert(L.begin(),L.end());
 }
};


class geo_nearest : public geowin_redraw, public geowin_update<std::list<Point>, std::list<Point> >
{
public:
  std::list<Point> pls,plt;

  virtual ~geo_nearest() {}
  
  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  { 
   std::list<Point>::const_iterator iter = pls.begin(), iter2 = plt.begin();
   for(;iter != pls.end();iter++,iter2++){
      leda_point pf = (*iter).to_float();
      leda_point pf2 = (*iter2).to_float();
      W.draw_arrow(pf.xcoord(),pf.ycoord(),pf2.xcoord(),pf2.ycoord(),c1);
   }
  }

  virtual bool write_postscript(ps_file& PS,leda_color c1,leda_color c2)
  {
   std::list<Point>::const_iterator iter = pls.begin(), iter2 = plt.begin();
   for(;iter != pls.end();iter++,iter2++){
      leda_point pf = (*iter).to_float();
      PS.draw_arrow(pf.xcoord(),pf.ycoord(),pf.xcoord(),pf.ycoord(),c1);
   }  
   return false;
  }
  
  virtual void update(const std::list<Point>& L, std::list<Point>&)
  { 
    pls.clear(); plt.clear();
    std::list<Point>::const_iterator iter = L.begin();
    
    for(;iter != L.end(); iter++){
      Vertex_handle vh = dt.nearest_vertex(*iter);
      pls.push_back(*iter);
      plt.push_back(vh->point());
    }
  }  
};


class geo_triangles : public geowin_redraw, public geowin_update<std::list<Point>, std::list<Point> >
{
public:
  std::list<Triangle> LT;
  geo_scene lines;

  virtual ~geo_triangles() {}

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {  
    std::list<Triangle>::const_iterator it   = LT.begin(), stop = LT.end();
    leda_color old = W.set_fill_color(leda_green);
    while( it != stop )
    {
       W << CGAL::convert_to_leda(*it);
       it++;
    }
    W.set_fill_color(old);
  }
  
  virtual bool write_postscript(ps_file& PS,leda_color c1,leda_color c2)
  {
    std::list<Triangle>::const_iterator it   = LT.begin(), stop = LT.end();
    leda_color old = PS.set_fill_color(leda_green);
    while( it != stop )
    {
       PS << CGAL::convert_to_leda(*it); 
       it++;
    }
    PS.set_fill_color(old);  
    return false;
  }  

  virtual void update(const std::list<Point>& L, std::list<Point>&)
  { 
    LT.clear();

    if (dt.dimension() != 2) return;

    GeoWin* gw = get_geowin(lines);

    std::list<Line> LST;
    gw->get_objects(lines,LST);

    std::list<Line>::const_iterator it;
    Line lakt;
    Point p1,p2; 

    for(it=LST.begin(); it != LST.end(); ++it) { 
       lakt= *it;
       p1=lakt.point1(); p2=lakt.point2();
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


class geo_triangles2 : public geowin_redraw, public geowin_update<std::list<Point>, std::list<Point> >
{
public:
  std::list<Triangle> LT;
  geo_scene locate_points;

  virtual ~geo_triangles2() {}

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {  
    std::list<Triangle>::const_iterator it   = LT.begin(), stop = LT.end();
    leda_color old = W.set_fill_color(leda_red);
    while( it != stop )
    {
       W << CGAL::convert_to_leda(*it);
       it++;
    }
    W.set_fill_color(old);
  }
  
  virtual bool write_postscript(ps_file& PS,leda_color c1,leda_color c2)
  {
    std::list<Triangle>::const_iterator it   = LT.begin(), stop = LT.end();
    leda_color old = PS.set_fill_color(leda_red);
    while( it != stop )
    {
       PS << CGAL::convert_to_leda(*it); 
       it++;
    }
    PS.set_fill_color(old);  
    return false;
  }  

  virtual void update(const std::list<Point>& L, std::list<Point>&)
  { 
    LT.clear();

    if (dt.dimension() != 2) return;

    GeoWin* gw = get_geowin(locate_points);   

    std::list<Point> LST;
    gw->get_objects(locate_points,LST);

    std::list<Point>::const_iterator it;
    Point lakt;

    for(it=LST.begin(); it != LST.end(); ++it) { 
       lakt= *it;
       Face_handle f = dt.locate(lakt);
       if (f != NULL && !dt.is_infinite(f)) LT.push_back(dt.triangle(f));
    }   
  }
};


int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("LEDA-rat_point"));
  geowin_init_default_type((std::list<Line>*)0, leda_string("LEDA-rat_line"));  

  std::list<Point>  L, LOC, NN;
  std::list<Line>   CGLL;

  GeoWin GW("CGAL - Delaunay triangulation demo");
  GW.add_help_text(leda_string("CGAL_delaunay_triang_2"));

  geo_scene my_scene= GW.new_scene(L); 
  // another input scene for the lines ...
  
  geo_scene line_scene= GW.new_scene(CGLL);
  GW.set_color(line_scene, leda_grey2);
  
  // another input scene for the points for location algor. ...
  geo_scene pointloc_scene= GW.new_scene(LOC);
  GW.set_color(pointloc_scene, leda_blue);
  
  geo_scene p_scene= GW.new_scene(NN);
  GW.set_color(p_scene, leda_red);  
  
  geo_delau delaunay_triang;
  geo_scene res2 = GW.new_scene(delaunay_triang,delaunay_triang, my_scene, "Delaunay Triangulation");
  GW.set_color(res2, leda_blue);
  GW.set_line_width(res2, 2);

  geo_triangles TRS;
  geo_scene sc1 = GW.new_scene( TRS, TRS, my_scene, "Triangles"); 
  GW.set_color(sc1, leda_blue);
  GW.set_visible(sc1,true);
  TRS.lines = line_scene;

  geo_triangles2 TPT;
  geo_scene res3 = GW.new_scene( TPT, TPT, my_scene, "Triangles2"); 
  GW.set_color(res3, leda_red);
  TPT.locate_points = pointloc_scene;
  
  geo_nearest NST;
  geo_scene res4 = GW.new_scene( NST, NST, p_scene, "Nearest vertex");   
  GW.set_color(res4, leda_blue2);
  
  GW.add_dependence(line_scene,sc1);
  GW.add_dependence(pointloc_scene,res3);
  GW.add_dependence(my_scene,res4);
  
  GW.set_all_visible(true);

  GW.edit(my_scene);
  return 0;  
}

#endif

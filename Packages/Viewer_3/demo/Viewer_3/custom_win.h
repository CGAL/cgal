//#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include "triang_2.h"




void close_c(Fl_Widget* w, void* v)
{

  Fl_Window* win = (Fl_Window*) v;
  delete win;
  
}




void tri_cb(Fl_Widget* w, void* v)
{
  CGAL::Viewer_3* win = (CGAL::Viewer_3*) v;
  typedef Triangulation_2::Point Point;
  Triangulation_2* tr = new Triangulation_2();  
  CGAL::Vector_2<rep_t> disp( 425.0, 425.0 );
  CGAL::Random_points_in_square_2<Point,CGAL::Creator_uniform_2<double,Point> >
    g ( 500.0 );
  for (int i=0; i<100; i++)
    tr->insert(*g++ + disp);
  CGAL::Drawable_triangulation_2<Triangulation_2>* dtr = new CGAL::Drawable_triangulation_2<Triangulation_2>(*tr,win->get_color(1), CGAL::ORANGE,CGAL::RAW, 2);
  win->add_drawable(dtr, win->get_group() + 1);
  win->display();
}


void del_cb(Fl_Widget* w, void* v)
{
  CGAL::Viewer_3* win = (CGAL::Viewer_3*) v;
  typedef Delaunay_2::Point Point;
  Delaunay_2* tr = new Delaunay_2();  
  CGAL::Vector_2<rep_t> disp( 425.0, 425.0 );
  CGAL::Random_points_in_square_2<Point,CGAL::Creator_uniform_2<double,Point> >
    g ( 500.0 );
  for (int i=0; i<100; i++)
    tr->insert(*g++ + disp);
  CGAL::Drawable_triangulation_2<Delaunay_2>* dtr = new CGAL::Drawable_triangulation_2<Delaunay_2>(*tr,win->get_color(1),  CGAL::ORANGE,CGAL::RAW, 2);
  win->add_drawable(dtr, win->get_group() + 1);
  win->display();
}


void vor_cb(Fl_Widget* w, void* v)
{
  CGAL::Viewer_3* win = (CGAL::Viewer_3*) v;
  typedef Delaunay_2::Point Point;
  Delaunay_2* tr = new Delaunay_2();  
  CGAL::Vector_2<rep_t> disp( 425.0, 425.0 );
  CGAL::Random_points_in_square_2<Point,CGAL::Creator_uniform_2<double,Point> >
    g ( 500.0 );
  for (int i=0; i<100; i++)
    tr->insert(*g++ + disp);
  CGAL:: Drawable_voronoi_2<Delaunay_2>* vo = new CGAL::Drawable_voronoi_2<Delaunay_2>(*tr,win->get_color(1), CGAL::RAW, 2);
  win->add_drawable(vo, win->get_group() + 1);
  win->display();
}


void vord_cb(Fl_Widget* w, void* v)
{
  CGAL::Viewer_3* win = (CGAL::Viewer_3*) v;
  typedef Delaunay_2::Point Point;
  Delaunay_2* tr = new Delaunay_2();  
  CGAL::Vector_2<rep_t> disp( 425.0, 425.0 );
  CGAL::Random_points_in_square_2<Point,CGAL::Creator_uniform_2<double,Point> >
    g ( 500.0 );
  for (int i=0; i<100; i++)
    tr->insert(*g++ + disp);
  CGAL::Drawable_voronoi_2<Delaunay_2>* vor = new CGAL::Drawable_voronoi_2<Delaunay_2>(*tr,CGAL::ORANGE, CGAL::RAW, 2);
  win->add_drawable(vor, win->get_group() + 1);
  CGAL::Drawable_triangulation_2<Delaunay_2>* dtr = new
    CGAL::Drawable_triangulation_2<Delaunay_2>(*tr,CGAL::BLUE,
					       CGAL::GREEN ,CGAL::RAW, 2);
  win->add_drawable(dtr, win->get_group() );

  win->display();
}




void my_panel(CGAL::GL_win* win, CGAL::Viewer_3* view)
{


 Fl_Window* flwin = new Fl_Window(100,500,"Demo win");

 Fl_Button* tri = new Fl_Button(5,5,80,25,"Triangulation");
 tri->callback(tri_cb,(void *) view);

 Fl_Button* del = new Fl_Button(5,35,80,25,"Delaunay");
 del->callback(del_cb,(void *) view);

 
 Fl_Button* vor = new Fl_Button(5,65,80,25,"Voronoi");
 vor->callback(vor_cb,(void *) view);

 Fl_Button* vord = new Fl_Button(5,95,80,25,"Vor. and Del.");
 vord->callback(vord_cb,(void *) view);


 Fl_Button* close = new Fl_Button(5,470,80,25,"Done");
 close->callback(close_c,(void *) flwin);
    flwin->end();
    flwin->show();

}



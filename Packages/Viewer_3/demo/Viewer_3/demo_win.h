#include <CGAL/point_generators_2.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
//#include "triang_2.h"

typedef CGAL::Cartesian<double> rep_t;
typedef CGAL::Triangulation_euclidean_traits_2< rep_t > Ttraits;
typedef CGAL::Triangulation_vertex_base_2<Ttraits>      Vertex_base ;
typedef CGAL::Triangulation_face_base_2<Ttraits>        Face_base ;
typedef CGAL::Triangulation_default_data_structure_2<Ttraits,Vertex_base,Face_base> TDS ;
typedef CGAL::Delaunay_triangulation_2<Ttraits, TDS> Delaunay_2;
typedef Delaunay_2::Point Point;
Delaunay_2* tr; 


void step(int val)
{
  double a=5;
  for (int i=0; i!=val*1000;i++)
    a= (a*a)/3 + (a*a)/3;
}




void close_c(Fl_Widget* w, void* v)
{

  Fl_Window* win = (Fl_Window*) v;
  delete win;
  
}


void del_cb(Fl_Widget* w, void* v)
{
  CGAL::Viewer_3* win = (CGAL::Viewer_3*) v;
  tr->clear();
  CGAL::Vector_2<rep_t> disp( 425.0, 425.0 );
  CGAL::Random_points_in_square_2<Point,CGAL::Creator_uniform_2<double,Point> >
    g ( 500.0 );
  for (int i=0; i<50; i++)
    tr->insert(*g++ + disp);
  CGAL::Drawable_triangulation_2<Delaunay_2>* dtr = new CGAL::Drawable_triangulation_2<Delaunay_2>(*tr,win->get_color(1),  CGAL::ORANGE,CGAL::RAW, 2);
  win->add_drawable(dtr, win->get_group()+1);
  win->display();
}


void parab_cb(Fl_Widget* w, void* v)
{
  CGAL::Viewer_3* win = (CGAL::Viewer_3*) v;
  tr->clear();
  
  float i= 0.00000000000001 ;
  while (i < 0.0000000001) {
    tr->insert(Point(i, i*i));
    i=i+0.00000000000001;
  } 
  cerr << "number of faces : " << tr->number_of_faces() << endl;
  CGAL::Drawable_triangulation_2<Delaunay_2>* dtr = new CGAL::Drawable_triangulation_2<Delaunay_2>(*tr,win->get_color(1),  CGAL::ORANGE,CGAL::RAW, 2);
  win->add_drawable(dtr, win->get_group()+1);
  win->display();
}




void iterator_cb(Fl_Widget* w, void* v)
{
  CGAL::Viewer_3* win = (CGAL::Viewer_3*) v;
  Delaunay_2::Face_iterator it=tr->finite_faces_begin();
  int gr = win->get_group()+1;
  int i=1;
  int count = 1;
  for (; it != tr->finite_faces_end() ;it++)
    {
      Delaunay_2::Triangle trg = tr->triangle(it);
      CGAL::Drawable_triangle_2<Delaunay_2::Triangle>*
	tri= new
	CGAL::Drawable_triangle_2<Delaunay_2::Triangle>(trg, CGAL::ORANGE, CGAL::DEEPBLUE);

      win->add_drawable(tri,gr);

      win->display();

      //      step(1000);
      //      sleep(1);
      win->remove_drawable(gr, 1);
      if (i==500) {
	cerr << "nombre de faces visitees : " << count << endl;
	i=0;
      }
      count ++;
      i++;
    }
  cerr << "delete group" << endl;
  win->delete_group(gr);
}

void line_face_cb(Fl_Widget* w, void* v)
{
  CGAL::Viewer_3* win = (CGAL::Viewer_3*) v;
  CGAL::Vector_2<rep_t> disp( 425.0, 425.0 );
  CGAL::Random_points_on_segment_2<Point,CGAL::Creator_uniform_2<double,Point> >
    g1 (Point(-100,10),Point(-100,900));
  CGAL::Random_points_on_segment_2<Point,CGAL::Creator_uniform_2<double,Point> >
    g2 (Point(900,10),Point(900,900));
  int gr = win->get_group()+1;
   Delaunay_2::Segment sg=Delaunay_2::Segment(*g1,*g2);
  CGAL::Drawable_segment_2<Delaunay_2::Segment>*
      seg= new
      CGAL::Drawable_segment_2<Delaunay_2::Segment>(sg, CGAL::PURPLE);
 win->add_drawable(seg,gr);
  Delaunay_2::Line_face_circulator lfc=tr->line_walk(*g1,*g2);
  Delaunay_2::Line_face_circulator lfc_o=lfc;

  do {
   if (!(tr->is_infinite(lfc))) {
    Delaunay_2::Triangle trg = tr->triangle(lfc);
    CGAL::Drawable_triangle_2<Delaunay_2::Triangle>*
      tri= new
      CGAL::Drawable_triangle_2<Delaunay_2::Triangle>(trg, CGAL::ORANGE, CGAL::DEEPBLUE);

      win->add_drawable(tri,gr);

      win->display();

      step(1000);
   }
   lfc++;
  }

  while (lfc != lfc_o);
 }


void circulator_cb(Fl_Widget* w, void* v)
{
  
  CGAL::Viewer_3* win = (CGAL::Viewer_3*) v;
  Delaunay_2::Vertex_iterator vit=tr->finite_vertices_begin();
  Delaunay_2::Face_circulator fct, fct_o;
  int gr = win->get_group()+1;
  for (; vit != tr->finite_vertices_end() ;vit++)
    {
      fct=tr->incident_faces(vit);
      fct_o=fct;
      do 
	{
          if (!(tr->is_infinite(fct))) {
            Delaunay_2::Triangle trg = tr->triangle(fct);
	    CGAL::Drawable_triangle_2<Delaunay_2::Triangle>*
	      tri= new
	      CGAL::Drawable_triangle_2<Delaunay_2::Triangle>(trg, CGAL::ORANGE, CGAL::DEEPBLUE);

	    win->add_drawable(tri,gr);

	    win->display();

	    step(1000);

	    win->remove_drawable(gr, 1);
	
	  }
	  fct++;
	}
      while (fct != fct_o);
    }
  win->delete_group(gr);
}
void demo_panel(CGAL::GL_win* win, CGAL::Viewer_3* view)
{


 Fl_Window* flwin = new Fl_Window(100,500,"Demo win");

 // Fl_Button* tri = new Fl_Button(5,5,80,25,"Triangulation");
 // tri->callback(tri_cb,(void *) view);

  Fl_Button* del = new Fl_Button(5,35,80,25,"Delaunay");
  del->callback(del_cb,(void *) view);

 
 Fl_Button* iter = new Fl_Button(5,65,80,25,"iterator");
 iter->callback(iterator_cb,(void *) view);

 Fl_Button* circu = new Fl_Button(5,95,80,25,"circulator");
 circu->callback(circulator_cb,(void *) view);

 Fl_Button* lfc = new Fl_Button(5,125,80,25,"Line face");
 lfc->callback(line_face_cb,(void *) view);

  Fl_Button* para = new Fl_Button(5,155,80,25,"Parabola");
  para->callback(parab_cb,(void *) view);


 Fl_Button* close = new Fl_Button(5,470,80,25,"Done");
 close->callback(close_c,(void *) flwin);
    flwin->end();
    flwin->show();

}

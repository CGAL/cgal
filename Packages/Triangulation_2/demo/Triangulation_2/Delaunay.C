#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) && !defined(CGAL_USE_CGAL_WINDOW)
int main(int argc, char* argv[])
{

  std::cout << "Sorry, this demo needs LEDA for visualisation. or CGAL_WINDOW";
  std::cout << std::endl;

  return 0;
}

#else
#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_drawing_mode CGAL::drawing_mode
#define leda_xor_mode     CGAL::xor_mode
#define leda_src_mode     CGAL::src_mode
#define leda_red          CGAL::red
#define leda_yellow       CGAL::yellow
#define leda_invisible    CGAL::invisible
#endif


#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Window_stream.h>

typedef double                              coord_type;
typedef CGAL::Cartesian<coord_type>         Gt;
typedef CGAL::Delaunay_triangulation_2<Gt>  Delaunay;
typedef CGAL::Window_stream                 Window_stream;

typedef Gt::Point_2                         Point;
typedef Delaunay::Face_handle               Face_handle;
typedef Delaunay::Vertex_handle             Vertex_handle;
typedef Delaunay::Edge                      Edge;

#include "parse.h"
#include "common_stuff.h"



void show_nearest_vertex(Delaunay &T, Window_stream &W)
{
    if (T.dimension()<1) return;
    std::cerr << "The vertex that is nearest to the cursor is highlighted"
	      << std::endl;
    std::cerr << "Click any button to continue" << std::endl;

    Vertex_handle nv = NULL;
    Vertex_handle v = NULL;

    Point p;
    Point q(coord_type(W.xmin()-1),coord_type(W.ymin()-1));

    leda_drawing_mode dm = W.set_mode(leda_xor_mode);
    while(1) {
        double x, y;
        int b = W.get_mouse(x,y);

        p = Point(coord_type(x),  coord_type(y));
        if(p != q){
	  v = T.nearest_vertex(p);
	  dm=W.set_mode(leda_xor_mode);
	  W << p<<v->point();
	  v = T.nearest_vertex(p); 
	  W << p << v->point();
	  W.set_mode(dm);
        }

        if( (nv != (CGAL_NULL_TYPE) NULL) && 
	    ( (b != NO_BUTTON) || ((p != q) && 
				   (v != nv) ) )){
            // unhighlight
            x = CGAL::to_double(nv->point().x());
            y = CGAL::to_double(nv->point().y());
            W.draw_node(x, y, leda_red);
        }
        if(b != NO_BUTTON){
            W.set_mode(dm);
            return;
        }
        if( (p != q) && (v != nv) ){
            nv = v;
            x = CGAL::to_double(nv->point().x());
            y = CGAL::to_double(nv->point().y());
            W.draw_node(x, y, leda_red);
            q = p;
        }
    }
}

void show_conflicts( Point& p, Delaunay& T, Window_stream &W )
{
  W.set_mode(leda_src_mode);
  W.clear();
  W << CGAL::BLUE << T;
  W << CGAL::RED << p;
  W.set_fill_color(leda_yellow);
  if (T.dimension() == 2) {
    std::list<Face_handle> conflict_faces;
    std::list<Edge>  hole_bd;
    T.get_conflicts_and_boundary(p, 
				 std::back_inserter(conflict_faces),
				 std::back_inserter(hole_bd));
    std::list<Face_handle>::iterator fit = conflict_faces.begin();
    std::list<Edge>::iterator eit = hole_bd.begin();
    for( ; fit != conflict_faces.end(); fit++)  {
      if ( ! T.is_infinite(*fit)) 
	W << CGAL::BLUE << T.triangle(*fit);
    }
    for( ; eit != hole_bd.end(); eit++)  {
      if ( ! T.is_infinite(*eit))
	W << CGAL::RED << T.segment(*eit);
    }
    W << CGAL::RED << p;
    any_button(W);
    T.star_hole(p, 
		hole_bd.begin(), 
		hole_bd.end(), 
		conflict_faces.begin(),
		conflict_faces.end());
    T.is_valid();
    //T.insert(p);
    W << CGAL::BLUE << T;
    any_button(W);
    W.clear();
    W.set_fill_color(leda_invisible);
    W << CGAL::BLUE << T;
  }
  else  {
    T.insert(p);
    W.clear(); W << CGAL::BLUE << T;
  }
  return;
}

void show_dual( Delaunay &T, Window_stream &W )
{
   std::cerr << "The dual of the triangulation is displayed" << std::endl;
   W << CGAL::RED;
   Delaunay::Face_iterator fit, fbegin=T.faces_begin(), fend=T.faces_end();
   for (fit=fbegin; fit != fend; ++fit)
     W << T.dual(fit);

   Delaunay::Edge_iterator eit, ebegin=T.edges_begin(), eend=T.edges_end();
   for (eit=ebegin; eit != eend; ++eit)
     {
       CGAL::Object o = T.dual(eit);
       Gt::Ray_2 r;
       Gt::Segment_2 s;
       Gt::Line_2 l;
       if (CGAL::assign(s,o)) W << s;
       if (CGAL::assign(r,o)) W << r;
       if (CGAL::assign(l,o)) W << l;
     }
    any_button(W);
   T.draw_dual(W);
   any_button(W);
}
 



template <class TRIANGULATION>
void window_input(TRIANGULATION &T,
		  Window_stream &W)
{
  Point p;
  Vertex_handle hv;

  while(1) {
    std::cerr << "Enter points with the left button" << std::endl;
    std::cerr << "conflict zone is shown before insertions" << std::endl;
    std::cerr << "Remove points with the middle button" << std::endl;
    std::cerr << "Right button terminates input of points" << std::endl;
    double x, y;
    int b = W.read_mouse(x,y);
    p = Point(coord_type(x),
	      coord_type(y));
    switch(b){
    case MOUSE_BUTTON(1) :
      show_conflicts(p, T , W);
      W.clear(); W << T;
      break; 
    case MOUSE_BUTTON(2) :
      hv = T.nearest_vertex(p);
      T.remove(hv);
      W.clear();
      W << T;
      break;
    case MOUSE_BUTTON(3) :  
      return;
    }
  }
}



int main(int argc, char* argv[])
{
  Options opt;
  parse(argc, argv, opt);

  Window_stream W(opt.winx, opt.winy); // physical window size
  W.init(opt.min, opt.max, opt.min);   // logical window size
  W << CGAL::BLUE;
  CGAL::cgalize( W);
  W.display();

  Delaunay T;
  file_input(T, W, opt);
  W << T;
  window_input(T, W);
  show_nearest_vertex(T,W);
  show_dual(T, W);
 }

#endif // CGAL_USE_LEDA || CGAL_USE_CGAL_WINDOW

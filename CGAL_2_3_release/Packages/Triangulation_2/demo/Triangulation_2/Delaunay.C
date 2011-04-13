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
#define leda_xor_mode CGAL::xor_mode
#define leda_src_mode CGAL::src_mode
#define leda_red      CGAL::red
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
    W << CGAL::RED << p;
   if (T.dimension() < 2){
    T.insert(p);
  }
  else {
    std::list<Face_handle> conflict_faces;
    std::list<Edge>  hole_bd;
    T.get_conflicts_and_boundary(p, 
				 std::back_inserter(conflict_faces),
				 std::back_inserter(hole_bd));
    std::list<Face_handle>::iterator fit = conflict_faces.begin();
    std::list<Edge>::iterator eit = hole_bd.begin();
    W << CGAL::RED ;
    for( ; fit != conflict_faces.end(); fit++)  {
      draw_face( *fit, T, W);
    }
    W << CGAL::BLACK;
    for( ; eit != hole_bd.end(); eit++)  {
      draw_edge( *eit, T, W);
    }
    any_button(W);
    T.insert(p);
    W << CGAL::BLUE << T;
    W << CGAL::BLACK;
    for(eit=hole_bd.begin() ; eit != hole_bd.end(); eit++)  {
      draw_edge( *eit, T, W);
    }
    any_button(W);
  }
  W.clear();
  W << CGAL::BLUE << T;
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
		  Window_stream &W,
		  const Options& opt)
{
    std::cerr << "Enter points with the left button" << std::endl;
    std::cerr << "conflict zone is shown before insertions" << std::endl;
    std::cerr << "Remove points with the middle button" << std::endl;
    std::cerr << "Right button terminates input of points" << std::endl;

    Point p;
    Point q(coord_type(W.xmin()-1),
            coord_type(W.ymin()-1));

    Face_handle highlight = NULL;
    Vertex_handle hv;

    while(1) {
        double x, y;
        int b = W.get_mouse(x,y);
        bool button_pressed = (b == MOUSE_BUTTON(1)) ||
                              (b == MOUSE_BUTTON(2)) ||
                              (b == MOUSE_BUTTON(3));
        p = Point(coord_type(x),
                  coord_type(y));
        bool mouse_moved = p != q;
        bool face_change = true,
             vertex_change = false;
        if( (highlight != (CGAL_NULL_TYPE) NULL) && 
	                  (button_pressed || mouse_moved) ){
            face_change = mouse_moved &&
                          ( T.oriented_side(highlight, p)
                                == CGAL::ON_NEGATIVE_SIDE );
            vertex_change = face_change ||
                            ( mouse_moved &&
                              ( hv != closest_vertex(T,highlight, p)));

            leda_drawing_mode dm = W.set_mode(leda_xor_mode);
            if(vertex_change){
                W << CGAL::RED ;
                W.draw_node(CGAL::to_double(hv->point().x()),
                            CGAL::to_double(hv->point().y()));
            }
            W.set_mode(leda_src_mode);
            if(face_change){
                W << CGAL::BLUE << T.triangle(highlight);
                highlight = NULL;
            }
            W.set_mode(dm);
        }
        if(b == MOUSE_BUTTON(1)){
            typename TRIANGULATION::Locate_type lt;
	    show_conflicts(p, T , W);
	    if(opt.check){
                T.is_valid();
            }

            if(lt != TRIANGULATION::VERTEX){
                W.clear();
                W << T;
            }
        } else if(b == MOUSE_BUTTON(2)){
            if(hv != (CGAL_NULL_TYPE) NULL){
                T.remove(hv);
                face_change = vertex_change = true;
                highlight = NULL;
                W.clear();
                W << T;
            }
        } else if(b == MOUSE_BUTTON(3)){
            // we are done. Nothing is highlighted
            break;
        }
        if( button_pressed || face_change){
            bool outside = highlight == (CGAL_NULL_TYPE) NULL;
            highlight = T.locate(p, highlight);
            if((highlight != (CGAL_NULL_TYPE) NULL) && 
	                     (! T.is_infinite(highlight)) &&
	                        T.dimension()==2){
                vertex_change = outside && true;
                leda_drawing_mode dm = W.set_mode(leda_src_mode);
                W << CGAL::BLUE << T.triangle(highlight) << CGAL::BLUE;
                W.set_mode(dm);
            } else {
                highlight = NULL;
            }
        }
        if(vertex_change){
            hv.clear();
        }
        if(button_pressed || vertex_change){
            if((highlight != (CGAL_NULL_TYPE) NULL) && 
	       (! T.is_infinite(highlight))){
                leda_drawing_mode dm = W.set_mode(leda_xor_mode);
                W << CGAL::RED;
                hv = closest_vertex(T, highlight, p);
                W.draw_node(CGAL::to_double(hv->point().x()),
                            CGAL::to_double(hv->point().y()));
                W << CGAL::BLUE;
                W.set_mode(dm);
            }
        }
        q = p;
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
  //W.set_fill_color(leda_red);
  window_input(T, W, opt);
  show_nearest_vertex(T,W);
  show_dual(T, W);
 }

#endif // CGAL_USE_LEDA || CGAL_USE_CGAL_WINDOW

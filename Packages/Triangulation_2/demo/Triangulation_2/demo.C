#include <CGAL/basic.h>
#include <fstream>
#include <list>


#if !defined(CGAL_USE_LEDA) && !defined(CGAL_USE_CGAL_WINDOW)
int main(int argc, char* argv[])
{
  std::cout << "Sorry, this demo needs LEDA or CGAL::WINDOW for visualisation";
  std::cout << std::endl;
  return 0;
}

#else
#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_drawing_mode CGAL::drawing_mode
#define leda_xor_mode CGAL::xor_mode
#define leda_src_mode CGAL::src_mode
#define leda_red      CGAL::red
#define leda_pink     CGAL::pink
#endif //CGAL_USE_CGAL_WINDOW

//to get shorter names
#define Cartesian Ca

#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/IO/Window_stream.h>

#include "parse.h"


typedef double coord_type;
typedef CGAL::Cartesian<coord_type>  Gt;
typedef CGAL::Triangulation_2<Gt>  Triangulation;


typedef Triangulation::Face  Face;
typedef Triangulation::Vertex Vertex;
typedef Triangulation::Edge   Edge;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Face_handle Face_handle;

typedef Triangulation::Face_circulator  Face_circulator;
typedef Triangulation::Vertex_circulator  Vertex_circulator;

typedef Triangulation::Locate_type Locate_type;

typedef Triangulation::Face_iterator  Face_iterator;
typedef Triangulation::Vertex_iterator  Vertex_iterator;
typedef Triangulation::Edge_iterator  Edge_iterator;
typedef Triangulation::Line_face_circulator  Line_face_circulator;
typedef Triangulation::Edge_circulator  Edge_circulator;
typedef Gt::Point_2          Point;
typedef Gt::Segment_2        Segment;
typedef Gt::Line_2           Line;
typedef CGAL::Window_stream  Window_stream;

#include "common_stuff.h"

template <class TRIANGULATION>
void window_input(TRIANGULATION &tr,
		  Window_stream &win)
{
  std::cerr << "Enter points with the left button" << std::endl;
  std::cerr << "Remove points with the middle button" << std::endl;
  std::cerr << "Right button terminates input of points" << std::endl;

  Point p;
  Vertex_handle vh;
  Face_handle  fh;

  while(1) {
    double x, y;
    int b = win.read_mouse(x,y);
    p = Point(coord_type(x),
	      coord_type(y));
    switch(b){
    case MOUSE_BUTTON(1) :
      tr.insert(p);
      win.clear();
      win << tr;
      break; 
    case MOUSE_BUTTON(2) :
      fh = tr.locate(p);
      vh = closest_vertex(tr, fh, p );
      tr.remove(vh);
      win.clear();
      win << tr;
      break;
    case MOUSE_BUTTON(3) :  
      return;
    }
  }
}


void draw_incident_edges(Triangulation &T,
                    Vertex_handle v,
                    Window_stream &W)
{
    leda_drawing_mode dm = W.set_mode(leda_xor_mode);

    W << CGAL::RED;
    Point p = v->point();

    Vertex_circulator vc = v->incident_vertices(),
                      done(vc);
    do {
        if(! T.is_infinite(vc)) {
            W << Segment(p, vc->point());
        }
    } while(++vc != done);

    W.set_mode(dm);
}

void redraw_incident_edges(Triangulation &T,
                    Vertex_handle v,
                    Window_stream &W)
{
    leda_drawing_mode dm = W.set_mode(leda_xor_mode);

    W << CGAL::RED;

    Edge_circulator ec = v->incident_edges(),
                      done(ec);
    do {
        if(! T.is_infinite(*ec)) {
            W << T.segment(*ec);
        }
    } while(++ec != done);

    W.set_mode(dm);
}

void draw_incident_edges(Triangulation &T,
                         Window_stream &W)
{
    if (T.dimension()<1) return;
    std::cerr << "Select a face" << std::endl;
    std::cerr << "The adjacent edges of the vertices of this face will be"
		 << std::endl
         << "highlighted in turn."<< std::endl;

    Point p;
    W >> p;
    Face_handle f = T.locate(p);

    for(int j = 0; j < 3; j++) {
        Vertex_handle v =  f->vertex(j);
        if(! T.is_infinite(v)) {
            std::cerr << "degree(v) = " << v->degree() << std::endl;
            draw_incident_edges(T, v, W);
            any_button(W);
            redraw_incident_edges(T, v, W);
        }
    }
}

void draw_faces_along_line(Triangulation &T,
                           Window_stream &W)
{
    if (T.dimension()<2) return;
    Point p, q;
    std::cerr << "Enter two points" << std::endl;
    std::cerr << "The faces intersected by the line joining those points "
	      << std::endl;
    std::cerr << " will be highlighted" << std::endl;
    std::cerr << std::endl;
    W << CGAL::RED;
    leda_drawing_mode dm = W.set_mode(leda_xor_mode);
    W >> p >> q;
    while (p==q) W << q;
    W << p << q << Line(p,q);
    W.set_mode(dm);

    Face_handle f = T.locate(p);
    Line_face_circulator lfc = T.line_walk(p, q, f),
                         done(lfc);
    if(lfc == (CGAL_NULL_TYPE) NULL){
        std::cerr << "Line does not intersect convex hull" << std::endl;
    } else {
      W.set_fill_color(leda_pink);
        do{
            if(! T.is_infinite( lfc  )){
                W << T.triangle( lfc );
            }
        }while(++lfc != done);
	//redraw the line
	W << p << q << Line(p,q);
	
        any_button(W);
	//clear all
	W.clear();
	W << CGAL::BLUE << T ;
    }
}

void draw_convex_hull(Triangulation &T,
                 Window_stream &W)
{
    if (T.dimension()<1) return;
    Point p, q;
    std::cerr << "Highlighting of the convex hull"<< std::endl;

    Vertex_circulator chc = T.infinite_vertex()->incident_vertices(),
                           done(chc);
    if(chc == (CGAL_NULL_TYPE) NULL) {
        std::cerr << "convex hull is empty" << std::endl;
    } else {
        leda_drawing_mode dm = W.set_mode(leda_src_mode);
        W << CGAL::RED;
        p = chc->point();
        do {
            --chc;
            q = chc->point();
            W << Segment(p, q);
            p = q;
        } while(chc != done);
        any_button(W);
        W << CGAL::BLUE;
        p = chc->point();
        do{
            ++chc;
            q = chc->point();
            W << Segment(p, q);
            p = q;
        } while(chc != done);
        W.set_mode(dm);
    }
}

void show_faces_iterator(Triangulation &T,
                 Window_stream &W)
{
  std::cerr << "Highlighting in turn each face traversed by the face iterarot "
			<<std::endl;
  W << CGAL::GREEN;
  Face_iterator fit= T.faces_begin();
  while (fit != T.faces_end()) {
    W << T.triangle(fit);
    fit++;
    any_button(W);
  }
  any_button(W);
  W << CGAL::BLUE;
  W << T;
  any_button(W);
}






void fileIO(Triangulation &T,
            Window_stream &W,
            const Options& opt)
{
    std::cerr << "The triangulation will be written to a file and read again"
			  << std::endl;
    {
        std::ofstream out("tr");
        CGAL::set_ascii_mode(out);
        out << T << std::endl;
    }
    Triangulation T2;

    std::ifstream in("tr");
    CGAL::set_ascii_mode(in);
    in >> T2;
    T2.is_valid();

    Window_stream W2(opt.winx, opt.winy);
    W2.init(opt.min, opt.max, opt.min);
    W2.set_show_coordinates(true);
    W2 << CGAL::BLUE;
    W2.display();

    W2 << T2;

    any_button(W);
}


int main(int argc, char* argv[])
{
    Options opt;
    parse(argc, argv, opt);
    Triangulation T, Copy;
  
    Window_stream W(opt.winx, opt.winy); // physical window size
    W.init(opt.min, opt.max, opt.min);   // logical window size
    W << CGAL::BLUE;
    CGAL::cgalize( W);
    W.display();
    
    file_input(T, W, opt);
    W << T;
    window_input(T, W);

    //container_input(T, W);
    draw_faces_along_line(T, W);
    draw_incident_edges(T, W);
    draw_convex_hull(T, W);
    // The show_faces_iterator is not so nice using any-button
    // and sleep is not portable 
    //show_faces_iterator(T,W);
 

    fileIO(T, W, opt);
    return 0;
}

#endif // CGAL_USE_LEDA || CGAL_USE_CGAL_WINDOW

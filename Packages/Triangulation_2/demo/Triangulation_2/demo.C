#include <CGAL/basic.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream>

//to get shorter names
#define Cartesian Ca

#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h>   // to avoid a g++ problem
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Segment_2.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/IO/Window_stream.h>

#include "parse.h"


typedef double coord_type;
typedef CGAL::Cartesian<coord_type>  Rp;

typedef CGAL::Point_2<Rp>  Point_;
typedef CGAL::Segment_2<Rp>  Segment_;
typedef CGAL::Ray_2<Rp>  Ray_;
typedef CGAL::Line_2<Rp>  Line_;
typedef CGAL::Triangle_2<Rp>  Triangle_;

typedef CGAL::Triangulation_euclidean_traits_2<Rp> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>  Triangulation_;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>  Delaunay_;

typedef Triangulation_::Face  Face_;
typedef Triangulation_::Vertex Vertex_;
typedef Vertex_::Vertex_handle Vertex_handle_;
typedef Face_::Face_handle Face_handle_;

typedef Triangulation_::Face_circulator  Face_circulator_;
typedef Triangulation_::Vertex_circulator  Vertex_circulator_;

typedef Triangulation_::Locate_type Locate_type_;

typedef Triangulation_::Face_iterator  Face_iterator_;
typedef Triangulation_::Vertex_iterator  Vertex_iterator_;
typedef Triangulation_::Edge_iterator  Edge_iterator_;
typedef Triangulation_::Line_face_circulator  Line_face_circulator_;
typedef Triangulation_::Edge_circulator  Edge_circulator_;

typedef CGAL::Window_stream  Window_stream;

#ifdef __GNU__
template < class R >
bool operator<(const CGAL::Point_2<R>& p, const CGAL::Point_2<R>& q)
{
    return CGAL::compare_lexicographically_xy (p, q) == CGAL::SMALLER;
}
#endif // __GNU__

Window_stream *W_global;
Triangulation_ *T_global;

void
any_button(CGAL::Window_stream &W)
{
    double x, y;
    std::cerr << "Press any button to continue" << std::endl;
    W.read_mouse(x,y);
}


template < class TRIANGULATION > 
Vertex_handle_ closest_vertex(const TRIANGULATION &T,
               Face_handle_ f,
               const Point_& p)
{
    Vertex_handle_ v = f->vertex(0);
    Rp::FT d  = CGAL::squared_distance(p, v->point());
    Rp::FT d2 = CGAL::squared_distance(p, f->vertex(1)->point());
    if(d2 < d){
        d = d2;
        v = f->vertex(1);
    }
    d2 = CGAL::squared_distance(p, f->vertex(2)->point());
    if(d2 < d){
        v = f->vertex(2);
    }
    return v;
}

template <class TRIANGULATION>
void window_input(TRIANGULATION &T,
		  Window_stream &W,
		  const Options& opt)
{
    std::cerr << "Enter points with the left button" << std::endl;
    std::cerr << "Remove points with the middle button" << std::endl;
    std::cerr << "Right button terminates input of points" << std::endl;

    Point_ p;
    Point_ q(coord_type(W.xmin()-1),
            coord_type(W.ymin()-1));

    Face_handle_ highlight = NULL;
    Vertex_handle_ hv;

    while(1) {
        double x, y;
        int b = W.get_mouse(x,y);
        bool button_pressed = (b == MOUSE_BUTTON(1)) ||
                              (b == MOUSE_BUTTON(2)) ||
                              (b == MOUSE_BUTTON(3));
        p = Point_(coord_type(x),
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
                              ( hv != closest_vertex(T, highlight, p)));

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
	    int li;
	    Face_handle_ loc = T.locate(p, lt, li);
            T.insert(p,lt, loc,li) ;
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
                W << CGAL::RED << T.triangle(highlight) << CGAL::BLUE;
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

void file_input(Triangulation_ &T,
           Window_stream& W,
           const Options& opt)
{
    if(! opt.file_input){
        return;
    }

    std::ifstream is(opt.fname);
    CGAL::set_ascii_mode(is);

    int n, count = 0;
    is >> n;
    std::cerr << "Reading " << n << " points" << std::endl;

    if( (! opt.check) && (! opt.draw)){
	std::cerr << "Reading from iterator" << std::endl;
        std::istream_iterator<Point_> begin(is);
        std::istream_iterator<Point_> end;
        T.insert(begin, end);
    }else{
	std::cerr << "Reading from point to point" << std::endl;
        Point_ mp;

        for(; n > 0; n--){
            is >> mp;
            T.insert(mp);
            if(opt.check){
		std::cerr << "Checking validity" << std::endl;
                T.is_valid();
            }
            if(opt.draw){
                W.clear();
                W << CGAL::BLUE << T << CGAL::RED;
            }

            if(++count == 100){
                std::cerr << ".";
                count = 0;
            }
        }
    }
    std::cerr << "Done with file input" << std::endl;
}

void container_input(Triangulation_ &T,
                Window_stream &W)
{
    std::list<Point_> L;
    L.push_front(Point_(0,0));
    L.push_front(Point_(1,0));
    L.push_front(Point_(1,1));

    int n = T.insert(L.begin(), L.end());
    std::cerr << n << " points inserted from a list." << std::endl;

    std::vector<Point_> V(3);
    V[0] = Point_(0, 0);
    V[1] = Point_(0.4, 0.4);
    V[2] = Point_(0.3, 0.3);

    n = T.insert(V.begin(), V.end());
    std::cerr << n << " points inserted from a vector." << std::endl;

    W.clear();
    W << T;
}

void draw_incident_edges(Triangulation_ &T,
                    Vertex_handle_ v,
                    Window_stream &W)
{
    leda_drawing_mode dm = W.set_mode(leda_xor_mode);

    W << CGAL::RED;
    Point_ p = v->point();

    Vertex_circulator_ vc = v->incident_vertices(),
                      done(vc);
    do {
        if(! T.is_infinite(vc)) {
            W << Segment_(p, vc->point());
        }
    } while(++vc != done);

    W.set_mode(dm);
}

void redraw_incident_edges(Triangulation_ &T,
                    Vertex_handle_ v,
                    Window_stream &W)
{
    leda_drawing_mode dm = W.set_mode(leda_xor_mode);

    W << CGAL::RED;

    Edge_circulator_ ec = v->incident_edges(),
                      done(ec);
    do {
        if(! T.is_infinite(*ec)) {
            W << T.segment(*ec);
        }
    } while(++ec != done);

    W.set_mode(dm);
}

void draw_incident_edges(Triangulation_ &T,
                         Window_stream &W)
{
    if (T.dimension()<1) return;
    std::cerr << "Select a face" << std::endl;
    std::cerr << "The adjacent edges of the vertices of this face will be"
		 << std::endl
         << "highlighted in turn."<< std::endl;

    Point_ p;
    W >> p;
    Face_handle_ f = T.locate(p);

    for(int j = 0; j < 3; j++) {
        Vertex_handle_ v =  f->vertex(j);
        if(! T.is_infinite(v)) {
            std::cerr << "degree(v) = " << v->degree() << std::endl;
            draw_incident_edges(T, v, W);
            any_button(W);
            redraw_incident_edges(T, v, W);
        }
    }
}

void draw_faces_along_line(Triangulation_ &T,
                           Window_stream &W)
{
    if (T.dimension()<2) return;
    Point_ p, q;
    std::cerr << "Enter two points" << std::endl;
    std::cerr << "The faces intersected by the line joining those points "
	      << std::endl;
    std::cerr << " will be highlighted" << std::endl;
    std::cerr << std::endl;
    W << CGAL::RED;
    leda_drawing_mode dm = W.set_mode(leda_xor_mode);
    W >> p >> q;
    while (p==q) W << q;
    W << p << q << Line_(p,q);
    W.set_mode(dm);

    Face_handle_ f = T.locate(p);
    Line_face_circulator_ lfc = T.line_walk(p, q, f),
                         done(lfc);
    if(lfc == (CGAL_NULL_TYPE) NULL){
        std::cerr << "Line does not intersect convex hull" << std::endl;
    } else {
        do{
            if(! T.is_infinite( lfc  )){
                W << T.triangle( lfc );
            }
        }while(++lfc != done);

        any_button(W);
        // Remove the line and unhighlight again
        dm = W.set_mode(leda_xor_mode);
        W << p << q << Line_(p,q);
	W.set_mode(dm);
	W << CGAL::BLUE;
	do{
            if(! T.is_infinite( lfc )){
                W << T.triangle( lfc );
            }
        }while(--lfc != done);
	
    }
}

void draw_convex_hull(Triangulation_ &T,
                 Window_stream &W)
{
    if (T.dimension()<1) return;
    Point_ p, q;
    std::cerr << "Highlighting of the convex hull"<< std::endl;

    Vertex_circulator_ chc = T.infinite_vertex()->incident_vertices(),
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
            W << Segment_(p, q);
            p = q;
        } while(chc != done);
        any_button(W);
        W << CGAL::BLUE;
        p = chc->point();
        do{
            ++chc;
            q = chc->point();
            W << Segment_(p, q);
            p = q;
        } while(chc != done);
        W.set_mode(dm);
    }
}

void show_faces_iterator(Triangulation_ &T,
                 Window_stream &W)
{
  std::cerr << "Highlighting in turn each face traversed by the face iterarot "
			<<std::endl;
  W << CGAL::GREEN;
  Face_iterator_ fit= T.faces_begin();
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





void show_nearest_vertex(Delaunay_ &T,
                    Window_stream &W)
{
    if (T.dimension()<1) return;
    std::cerr << "The vertex that is nearest to the cursor is highlighted"
	          << std::endl;
    std::cerr << "Click any button to continue" << std::endl;

    Vertex_handle_ nv = NULL;
    Vertex_handle_ v = NULL;

    Point_ p;
    Point_ q(coord_type(W.xmin()-1),
            coord_type(W.ymin()-1));

    leda_drawing_mode dm = W.set_mode(leda_xor_mode);
    while(1) {
        double x, y;
        int b = W.get_mouse(x,y);

        p = Point_(coord_type(x),
                  coord_type(y));
        if(p != q){
	  v = T.nearest_vertex(p);dm=W.set_mode(leda_xor_mode);
	  W << p<<v->point();
	  v = T.nearest_vertex(p); W << p<<v->point();W.set_mode(dm);
        }

        if( (nv != (CGAL_NULL_TYPE) NULL) && ( (b != NO_BUTTON) || ((p != q) && (v != nv) ) )){
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

void fileIO(Delaunay_ &T,
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
    Triangulation_ T2;

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


void draw_dual( Delaunay_ &T, Window_stream &W )
{
   std::cerr << "The dual of the triangulation is displayed" << std::endl;
   W << CGAL::RED;
    Delaunay_::Face_iterator fit, fbegin=T.faces_begin(), fend=T.faces_end();
    for (fit=fbegin; fit != fend; ++fit)
        W << T.dual(fit);

    Delaunay_::Edge_iterator eit, ebegin=T.edges_begin(), eend=T.edges_end();
    for (eit=ebegin; eit != eend; ++eit)
    {
        CGAL::Object o = T.dual(eit);
        Gt::Ray r;
        Gt::Segment s;
        if (CGAL::assign(s,o)) W << s;
        if (CGAL::assign(r,o)) W << r;
    }
    any_button(W);
}
 
void action(int n)
{
  std::cerr << "action is called "<<n<<std::endl;
}


int main(int argc, char* argv[])
{
    Options opt;
    parse(argc, argv, opt);
    Triangulation_ T, Copy;
    Delaunay_ D, DCopy;
    T_global = &T;

    Window_stream W(opt.winx, opt.winy); // physical window size
    W_global = &W;

    W.init(opt.min, opt.max, opt.min);   // logical window size
    //  W.set_show_coordinates(true);
    W << CGAL::BLUE;
    //W.set_mode(leda_src_mode);
    //W.set_node_width(3);
    // W.button("toto",1,action);
    CGAL::cgalize( W);
    W.display();
    
     file_input(T, W, opt);
    (void) T.is_valid();
    W << T;
    window_input(T, W, opt);
    W.set_mode(leda_src_mode);
    W << CGAL::RED;
    W.clear();
    if(opt.check) T.is_valid();
    W << T;
    any_button(W);
    // Copy constructor
    std::cout << "copy"<<std::endl;
    Copy = Triangulation_(T);
    W << CGAL::BLUE;
    W << Copy;
    container_input(T, W);
    draw_faces_along_line(T, W);
    draw_incident_edges(T, W);
    draw_convex_hull(T, W);
    show_faces_iterator(T,W);
 
    std::cout <<std::endl<<std::endl<< "DELAUNAY TRIANGULATION"<<std::endl;
    T_global = (Triangulation_ *)&D;
    W.clear();
    window_input(D, W, opt);
    W.set_mode(leda_src_mode);
    W << CGAL::RED;
    W.clear();
    if(opt.check) D.is_valid();
    W << D;
    any_button(W);
    // Copy constructor
    std::cout << "copy"<<std::endl;
    DCopy = Delaunay_(D);
    W << CGAL::BLUE;
    W << DCopy;
    show_nearest_vertex(D, W);
    fileIO(D, W, opt);
    draw_dual(D, W);

    return 0;
}


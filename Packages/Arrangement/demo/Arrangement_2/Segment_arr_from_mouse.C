// demo/Arrangement_2/Segment_arr_from_mouse.C
//
//constructs a segment arrangement from CGAL window.
// We use the leda traits (therefore we are using leda functions).

#include <CGAL/config.h> // needed for the LONGNAME flag

#ifdef CGAL_CFG_NO_LONGNAME_PROBLEM
// Define shorter names to please linker (g++/egcs)
#define Arrangement_2 Ar
#define Arr_leda_segment_exact_traits Alset
#define Arr_2_default_dcel A2d
#define In_place_list_iterator IPLI
#define Arr_2_vertex_base Avb
#define Arr_2_halfedge_base Ahb
#define Arr_2_face_base Afb
#define Point_2 pT
#define Segment_2 sT
#define Topological_map TpM
#define _List_iterator Lit
#define Halfedge hE
#define Forward_circulator_tag Fct
#endif

#include <CGAL/basic.h>

#ifndef CGAL_USE_LEDA
int main()
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/Arr_leda_segment_exact_traits.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>

#ifndef CGAL_IO_FILE_DRAWER_H
#include <CGAL/IO/Pm_drawer.h>
#endif

#ifndef CGAL_IO_DRAW_PM_H
#include <CGAL/IO/draw_pm.h>
#endif

#include <CGAL/IO/Window_stream.h>

#include <vector>

typedef CGAL::Arr_leda_segment_exact_traits         Traits;

typedef Traits::Point                                  Point;
typedef Traits::X_curve                                X_curve;

typedef CGAL::Arr_base_node<X_curve>   Base_node;
typedef CGAL::Arr_2_default_dcel<Traits> Dcel;

typedef CGAL::Arrangement_2<Dcel,Traits,Base_node > Arr_2;

//I had to add these in global namespace for the program to compile
CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
                          const Point& p)
{
  //return os << leda_point(p.xcoordD(),p.ycoordD());
  return os << p.to_point();
}

CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
                          const X_curve &c)
{
  return os << c.to_segment();
}

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static Arr_2 arr;
static CGAL::Window_stream W(400, 400, "CGAL - Segment Arrangement Demo");


CGAL_BEGIN_NAMESPACE

class My_Arr_drawer : public Pm_drawer<Arr_2,Window_stream> {
private:
  typedef Pm_drawer<Arr_2,Window_stream>  Base;
public:
  My_Arr_drawer( Window_stream& W ): Pm_drawer<Arr_2,Window_stream>( W ){}
  
  void draw_face(Face_handle f) {
    if (f->does_outer_ccb_exist()) {
      Arr_2::Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
	W << cc->curve();
      } while (++cc != f->outer_ccb());  
    }

    Arr_2::Holes_iterator hit=f->holes_begin(),eit=f->holes_end();
    for (;hit!=eit; ++hit) {
      Arr_2::Ccb_halfedge_circulator cc=*hit; 
      do {
	W << cc->curve();
	} while (++cc != *hit);  
    }      
  }

  void draw_vertices(Vertex_const_iterator Vertices_begin, 
		     Vertex_const_iterator Vertices_end) {
    W << GREEN;
    Base::draw_vertices(Vertices_begin, Vertices_end);
  }
  
  void draw_halfedges(Halfedge_const_iterator Halfedges_begin, 
		      Halfedge_const_iterator Halfedges_end) {
    W << BLUE;
    Base base(window());
    base.draw_halfedges(Halfedges_begin, Halfedges_end);
  }
  
};
 
Window_stream& operator<<(Window_stream& os, Arr_2 &A)
{
  My_Arr_drawer drawer(os);
  
  draw_pm(arr, drawer, os);
  
  return os;
}

CGAL_END_NAMESPACE

// redraw function for the LEDA window. 
// used automatically when window reappears.
void redraw(CGAL::Window_stream * wp) 
{ wp->start_buffering();
  wp->clear();
  // draw arragnement
  *wp << arr;
  wp->flush_buffer();
  wp->stop_buffering();
}


int main()
{
  double x0=-400,x1=400,y0=-400;

  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.button("finish",10);
  W.open_status_window();
  W.display();

  //read input from window
  std::cout << "Left button to start and end the segment.\n";
  std::cout << "Clicking close to a vertex assumes the location" 
	    << "is at the vertex"
	    << std::endl;

  std::vector<Point> cv1;

  Point pnt;
  bool begin=true;

  for (;;) {
    double x, y;
    int b = W.get_mouse(x,y);
    if (b==10) break;
    pnt = Point(x,y);
  
    if (b == MOUSE_BUTTON(1))
      {
        
        for(Arr_2::Vertex_iterator vi = arr.vertices_begin();
            vi != arr.vertices_end(); ++vi) {
          //we are using the leda sqr_dist func
          if ( pnt.sqr_dist(vi->point()) < ((x1-x0)/50)*((x1-x0)/50) )
            pnt=vi->point();
        }
        
        cv1.push_back(pnt);
        W << CGAL::BLACK;
        W << pnt;
        W << CGAL::GREEN;
        
        if (!begin) {
          if ( cv1[0] == cv1[1] ){
            //Error. Segment has a zero length.
             W.set_status_string("Error. Segment has a zero length.");
            redraw( &W );
          }
          else{ 
            arr.insert(X_curve(cv1[0],cv1[1]));
            W.set_status_string("  ");
            W << arr;
	  }
          cv1.clear();
        }
        begin=!begin;
      }
  }

  
  W << arr;
   
  // Point Location Queries
  std::cout << "\nEnter a point with left button." << std::endl;
  W << CGAL::RED;

  Point p;

  Arr_2::Halfedge_handle e;
  
  // if map is empty
  if (arr.halfedges_begin() == arr.halfedges_end()) 
    {
      std::cout << std::endl;
      std::cout << "No edges were inserted. Planar map is empty. Exiting.";
      std::cout << std::endl;
    }
  else {
    // map is not empty
    
    CGAL::My_Arr_drawer  drawer(W);
    for (; ; ) {
      
      double x,y;
      int b=W.read_mouse(x,y);
      if (b==10) break;
      else
	p=Point(x,y);

      W << arr;
    
      Arr_2::Locate_type lt;
      e = arr.locate(p,lt);
      
      //color the face on the screen
      Arr_2::Face_handle f=e->face();
      drawer.draw_face(f);
      
      /*if (f->does_outer_ccb_exist()) {
	Arr_2::Ccb_halfedge_circulator cc=f->outer_ccb();
	do {
	W << cc->curve();
	} while (++cc != f->outer_ccb());
	
	}
      
	Arr_2::Holes_iterator hit=f->holes_begin(),eit=f->holes_end();
	for (;hit!=eit; ++hit) {
	Arr_2::Ccb_halfedge_circulator cc=*hit; 
	do {
	W << cc->curve();
	} while (++cc != *hit);
      
	} */     
    }
  } // else of 'if map is empty'

  return 0;  
}

#endif // CGAL_USE_LEDA

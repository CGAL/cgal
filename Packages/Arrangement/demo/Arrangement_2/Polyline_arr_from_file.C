#include <CGAL/config.h> // needed for the LONGNAME flag

#ifdef CGAL_CFG_NO_LONGNAME_PROBLEM
// Define shorter names to please linker (g++/egcs)
#define Arrangement_2 Ar
#define Cartesian cR
#define Arr_polyline_traits ARPT
#define Arr_2_default_dcel A2d
#define In_place_list_iterator IPLI
#define Arr_2_vertex_base Avb
#define Arr_2_halfedge_base Ahb
#define Arr_2_face_base Afb
#define Td_X_trapezoid TdXt
#define Td_traits Tdt
#define Planar_map_traits_wrap Pmtw
#define PL_X_curve_plus PXcp
#define Point_2 pT
#define allocator aR
#define Arr_base_node Abn
#define Topological_map TpM
#define _Arr_face_circ Afc
#define _Pm_Halfedge PmH
#define _List_iterator Lit
#define Halfedge hE
#define Forward_circulator_tag Fct
#endif

//constructs a ployline arrangement from a file.
// We use the leda traits (therefore we use leda functions).

// PARTIALLY CHANGED
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_polyline_traits.h>

#ifndef CGAL_USE_LEDA
int main(int argc, char* argv[])
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/leda_real.h>
#include <CGAL/IO/Window_stream.h>

typedef leda_real                            NT;
typedef CGAL::Cartesian<NT>                  Rep;

typedef CGAL::Arr_polyline_traits<Rep>       Traits;

typedef Traits::Point                        Point;
typedef Traits::X_curve                      X_curve;
typedef Traits::Curve                        Curve;

typedef CGAL::Arr_base_node<X_curve>         Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>     Dcel;

typedef CGAL::Arrangement_2<Dcel,Traits,Base_node > Arr_2;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static Arr_2               arr; 
static CGAL::Window_stream W(400, 400, "CGAL - Polyline Arrangement Demo");

leda_point to_leda_pnt(Point p)
{
  return leda_point(p.x().to_double(), p.y().to_double());
}
CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
                          const  Point& p)                         
{
  // conversion to leda_point in order to show it on screen
  return os << to_leda_pnt(p);
}

// draw a polyline, with points as 'x's
CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
				const X_curve &c)
{
  X_curve::const_iterator sit, tit;
  sit = c.begin();
  tit = sit; tit++;
  W.draw_point(to_leda_pnt(*sit), leda_green); // draw first point
  for (; tit != c.end(); tit++, sit++)
    {
    // conversion to screen drawble segment
    os << leda_segment(to_leda_pnt(*sit), to_leda_pnt(*tit));
    W.draw_point(to_leda_pnt(*tit), leda_green);
    }
    
  return os;
}

// draw an arrangement with polyline traits
CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
                                Arr_2 &A)
{
   Arr_2::Halfedge_iterator it = A.halfedges_begin();
   
   os << CGAL::BLUE;
   while(it != A.halfedges_end()){

      os << (*it).curve();
      ++it; ++it;
    }
   it = A.halfedges_begin();

    os.set_flush( 1 );
    os.flush();

    return os;
}

void read_arr(Arr_2 & arr, char * filename)
{  
  std::ifstream file(filename);

  int num_polylines, num_x_curves ;
  double x,y; // only used to read file not used in computations
  X_curve polyline;

  file >> num_polylines;
  while (num_polylines--) {
    file >> num_x_curves;
    while (num_x_curves--) {
      file >> x >> y;
      Point s(x,y);
      polyline.push_back(s);
    }

    arr.insert(polyline);

    polyline.clear();
  }
}

// redraw function for the LEDA window. used automatically when window reappears
void redraw(CGAL::Window_stream * wp) 
{ wp->start_buffering();
  wp->clear();
  // draw arragnement
  *wp << arr;
  wp->flush_buffer();
  wp->stop_buffering();
}

int main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cout << "usage: Polyline_arr_from_file filename\n";
    exit(1);
  }

  double x0=-20,x1=300,y0=-20;
  enum { THE_BUTTON = 10 };
  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.button("  Quit  ", THE_BUTTON);
  W.open_status_window();
  W.display(leda_window::center,leda_window::center);

  // read arrangement from file
  read_arr(arr, argv[1]);
  W << arr;

 // (3) Point Location part
  Arr_2::Halfedge_handle e;
  double x,y; // only used in drawing, not computations
  Point pnt;

  W.set_status_string("Enter a point with left button. Press Quit to quit.");  
  W.set_button_label(THE_BUTTON, "  Quit  ");
  while (W.read_mouse(x,y) != THE_BUTTON) {
    pnt = Point(x,y);
    W << arr;
    
    Arr_2::Locate_type lt;
    e = arr.locate(pnt ,lt);

    //color the face on the screen
    W << CGAL::GREEN;
    Arr_2::Face_handle f=e->face();
    if (f->does_outer_ccb_exist()) {
      Arr_2::Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
	W << cc->curve();
      } while (++cc != f->outer_ccb());
    }
    for (Arr_2::Holes_iterator ho=f->holes_begin(),hoe=f->holes_end();
	 ho!=hoe; ++ho) {
      Arr_2::Ccb_halfedge_circulator cc=*ho; 
      do {
	W << cc->curve();
      } while (++cc != *ho);
    }
  } // location

  return 0;  
}

#endif // CGAL_USE_LEDA

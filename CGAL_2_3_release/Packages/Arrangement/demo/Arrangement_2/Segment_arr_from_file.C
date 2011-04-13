// demo/Arrangement_2/Segment_arr_from_file.C
//
//constructs a segment arrangement from file.
// We use a homogeneous representation with the cached traits.

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
#define Homogeneous Ho
#define Quotient Qu
#define _In_place_list_iterator IPLI
#define Arrangement_2 Ar
#define Arr_segment_exact_traits ASET
#define Arr_base_node ABN
#endif

//File format is (coordinates should be between [-400,400] to be in window):
//#number_of_segments
//#x1 y1 x2 y2
// ....
 


#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
//#include <CGAL/Cartesian.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Timer.h>

#include <vector>
#include <fstream>

#ifndef CGAL_USE_LEDA
int main()
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/leda_integer.h>
//#include <CGAL/leda_real.h>
//#include <CGAL/leda_rational.h>
//#include <CGAL/Double_eps.h>

#include <CGAL/Arr_2_default_dcel.h>
//#include <CGAL/Arr_segment_exact_cached_traits.h>
#include <CGAL/Arr_segment_exact_traits.h>

#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Draw_preferences.h>

typedef leda_integer NT;
//typedef double           NT;
//typedef leda_real        NT;
//typedef leda_rational    NT;
//typedef CGAL::Double_eps NT;

typedef CGAL::Homogeneous<NT>                      Rep;
//typedef CGAL::Cartesian<NT>                      Rep;

//typedef CGAL::Arr_segment_exact_cached_traits<Rep>   Traits;
typedef CGAL::Arr_segment_exact_traits<Rep>          Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;


typedef CGAL::Arr_base_node<X_curve>   Base_node;
typedef CGAL::Arr_2_default_dcel<Traits> Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node > Arr_2;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static Arr_2               arr; 
static CGAL::Window_stream W(400, 400, "CGAL - Segment Arrangement Demo");

CGAL_BEGIN_NAMESPACE
Window_stream& operator<<(Window_stream& os, Arr_2 &A)
{
  My_Arr_drawer< Arr_2,
                 Arr_2::Ccb_halfedge_circulator, 
                 Arr_2::Holes_iterator> drawer(os);
  
  draw_pm(arr, drawer, os);
  
  return os;
}
CGAL_END_NAMESPACE

void color_face(CGAL::Window_stream& W, Arr_2::Halfedge_handle e,
                CGAL::Color color) {
  //color the face on the screen
  W << color;
  Arr_2::Face_handle f=e->face();
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
// redraw function for the LEDA window. used automatically when window 
// reappears
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
  //for Double_eps
  //CGAL::set_eps(0.0001); 

  if (argc != 2) {
    std::cout << "usage: Segment_arr_from_file filename\n";
    exit(1);
  }

  CGAL::Timer insrt_t;

  double x0=-200,x1=200,y0=-200;

  W.init(x0,x1,y0);
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.button("finish",10);
  W.open_status_window();
  W.display();

  std::ifstream file(argv[1]);
  int num_curves;
  file >> num_curves;
  while (num_curves--) {
    //NT x,y;
    double x,y; //(actually int)
    file >> x >> y;
    NT nx(x),ny(y);
    Point s(nx,ny);
    file >> x >> y;
    NT mx(x),my(y);
    Point t(mx,my);
    
    insrt_t.start();
    arr.insert(X_curve(s,t));
    insrt_t.stop();
  }
  
  W << arr;
  
  //POINT LOCATION
  W.set_status_string( "Enter a point with left button." );
  Point p;
  
  Arr_2::Halfedge_handle e=arr.halfedges_begin();
  
  for (; ; ) {
    
    //instead of before
    double x,y;
    int b=W.read_mouse(x,y);
    if (b==10) break;
    else
      p=Point(x,y);

    color_face(W,e,CGAL::BLUE);


    //W << arr;
    
    Arr_2::Locate_type lt;
    e = arr.locate(p,lt);

    //color the face on the screen
    color_face(W,e,CGAL::GREEN);
  }

  return 0;  
}

#endif // CGAL_USE_LEDA

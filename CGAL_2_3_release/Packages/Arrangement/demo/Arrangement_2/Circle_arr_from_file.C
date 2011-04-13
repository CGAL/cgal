// demo/Arrangement_2/Circle_arr_from_file.C
//
//constructs an arrangement of circles from file

//File format is:
//#number_of_circles
//#x_center y_center squared_radius
//# ....

#include <CGAL/Cartesian.h>
#include <fstream>
#include <CGAL/Timer.h>
#include <CGAL/Arr_circles_real_traits.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>

#ifndef CGAL_USE_LEDA
int main()
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/leda_real.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Arr_circle_traits_Window_stream.h>
#include <CGAL/Draw_preferences.h>

typedef CGAL::Arr_circles_real_traits<leda_real> Traits; 

typedef Traits::Point                                  Point;
typedef Traits::Curve                                  Curve; 
typedef Traits::X_curve                                X_curve;

typedef CGAL::Arr_base_node<Curve>                     Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>               Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >    Arr_2;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
static Arr_2 arr;
static CGAL::Window_stream W(400, 400, "CGAL - Circle Arrangement Demo");

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
  CGAL::Timer insrt_t;

  int circles_num;
  if (argc<2) {
    std::cerr << "usage: Circle_arr_from_file <filename>" << std::endl;
    exit(1);
  }

  //read circles from file
  std::ifstream f(argv[1]);
  f >> circles_num;

  std::vector<Curve> circles;
  
  double max_r2=1,max_x=1,min_x=-1,min_y=-1; //for adjusting the window size
  while (circles_num--) {
    leda_real x,y,r2;
    f >> x >> y >> r2;
    circles.push_back(Curve(x,y,r2));

    double dx=CGAL::to_double(x);
    double dy=CGAL::to_double(y);
    double dr2=CGAL::to_double(r2);
    if (dr2 > max_r2) max_r2=dr2;
    if (dx > max_x) max_x=dx;
    if (dx < min_x) min_x=dx;
    if (dy < min_y) min_y=dy;
  }
  f.close();

  W.init(min_x-2*sqrt(max_r2),max_x+2*sqrt(max_r2),min_y-2*sqrt(max_r2));
  W.set_redraw(redraw);
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.button("finish",10);
  W.open_status_window();
  W.display();

  for (unsigned int i=0; i<circles.size(); ++i) {
    std::cout << "inserting circle " << i+1 << std::endl;
    insrt_t.start();
    arr.insert(circles[i]);
    insrt_t.stop();
  }

  W << arr;
   std::cout << "Total insertion time: " <<  insrt_t.time() << std::endl;

  //POINT LOCATION
  W.set_status_string( "Enter a query point with left mouse button" );
  Point p;

  Arr_2::Halfedge_handle e;
  
  for (;;) {
    double x,y;
    int b=W.read_mouse(x,y);
    if (b==10) break;
    else
      p=Point(x,y);

    W << arr;
    
    Arr_2::Locate_type lt;
    e = arr.locate(p,lt);

    Arr_2::Face_handle fh=e->face();
    //Arr_2::Ccb_halfedge_circulator cc(e);
    Arr_2::Ccb_halfedge_circulator cc;

    if (fh != arr.unbounded_face()) {
      cc=fh->halfedge_on_outer_ccb();
      do {
        W << cc->curve();
      } while (++cc != fh->halfedge_on_outer_ccb());
    }

    Arr_2::Holes_iterator hit=fh->holes_begin(), eit=fh->holes_end();
    for (;hit!=eit;++hit) {
      cc=*hit;
      do {
        W << cc->curve();
      } while (++cc != *hit);
    }
      
  }

  return 0;  
}

#endif // CGAL_USE_LEDA

// demo/Arrangement_2/Circle_arr_from_file.C
//
// constructs an arrangement of circles from file
// 
// File format is:
// #number_of_circles
// #x_center y_center squared_radius
// # ....

#include <CGAL/basic.h>

#ifndef CGAL_USE_LEDA
#include <iostream>
int main()
{
  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;
  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Timer.h>
#include <CORE/BigInt.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Conic_arc_2_Window_stream.h>
#include <CGAL/Draw_preferences.h>
#include <fstream>

typedef CORE::BigInt                            CfNT;
typedef CGAL::Cartesian<CfNT>                   Int_kernel;
typedef Int_kernel::Point_2                     Int_point_2;
typedef Int_kernel::Circle_2                    Int_circle_2;

typedef CORE::Expr                              CoNT;
typedef CGAL::Cartesian<CoNT>                   Alg_kernel;

typedef CGAL::Arr_conic_traits_2<Int_kernel,
                                 Alg_kernel>    Traits_2;

typedef Traits_2::Point_2                       Point_2;
typedef Traits_2::Curve_2                       Curve_2;
typedef Traits_2::X_monotone_curve_2            X_monotone_curve_2;

typedef CGAL::Arr_2_default_dcel<Traits_2>      Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits_2>      Arr_2;

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
{ 
  wp->start_buffering();
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

  std::vector<Int_circle_2> circles;
  
  double max_r2=1,max_x=1,min_x=-1,min_y=-1; //for adjusting the window size
  while (circles_num--) 
  {
    int   x,y,r2;
    f >> x >> y >> r2;
    
    Int_point_2  center = Int_point_2 (CfNT(x), CfNT(y));
    circles.push_back(Int_circle_2 (center, r2));

    double dx = x;
    double dy = y;
    double dr2 = r2;
    if (dr2 > max_r2) max_r2=dr2;
    if (dx > max_x) max_x=dx;
    if (dx < min_x) min_x=dx;
    if (dy < min_y) min_y=dy;
  }
  f.close();

  W.init(min_x-2*sqrt(max_r2),max_x+2*sqrt(max_r2),min_y-2*sqrt(max_r2));
  W.set_redraw(redraw);
  W.set_mode(CGAL_LEDA_SCOPE::src_mode);
  W.set_node_width(3);
  W.button("finish",10);
  W.open_status_window();
  W.display();

  for (unsigned int i=0; i<circles.size(); ++i) 
  {
    std::cout << "inserting circle " << i+1 << std::endl;
    insrt_t.start();
    arr.insert(Curve_2(circles[i]));
    insrt_t.stop();
  }

  W << arr;
   std::cout << "Total insertion time: " <<  insrt_t.time() << std::endl;

  //POINT LOCATION
  W.set_status_string( "Enter a query point with left mouse button" );
  Point_2 p;

  Arr_2::Halfedge_handle e;
  
  for (;;) 
  {
    double x,y;
    int b=W.read_mouse(x,y);
    if (b==10) break;
    else
      p=Point_2(x,y);

    W << arr;
    
    Arr_2::Locate_type lt;
    e = arr.locate(p,lt);

    Arr_2::Face_handle fh=e->face();
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

#endif

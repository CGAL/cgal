#include <CGAL/Cartesian.h>

#include <iostream>

#include <CGAL/Fixed_precision_nt.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Window_stream.h>

#ifdef CGAL_USE_CGAL_WINDOW
#define leda_src_mode CGAL::src_mode
#endif

typedef CGAL::Fixed_precision_nt coord_type;
static bool Fixed_precision_nt_init_result 
                =  CGAL::Fixed_precision_nt::init(2000.0);

typedef CGAL::Cartesian<coord_type>               Repclass;
typedef Repclass::Point_2                         Point_;
typedef CGAL::Delaunay_triangulation_2<Repclass>  Delaunay_;

int main()
{
    CGAL::force_ieee_double_precision();

    Delaunay_ D;
    CGAL::Window_stream W(200,200); // physical window size

    W.init(-1,1,-1);   // logical window size
    W << CGAL::BLUE;
    W.set_mode(leda_src_mode);
    W.set_node_width(3);
    W.display();
    
    std::cout << std::endl << std::endl
              << "DELAUNAY TRIANGULATION" << std::endl;
    
    while(1) {
      double x, y;
      int b = W.get_mouse(x,y);
      if ( ! ( (b == MOUSE_BUTTON(1)) ||
	       (b == MOUSE_BUTTON(2)) ||
	       (b == MOUSE_BUTTON(3))))
        continue;
      if (b != MOUSE_BUTTON(1))
        break;
      D.insert( Point_(coord_type(x), coord_type(y)));
      W.clear();
      W << D;
    }
    return 0;
}

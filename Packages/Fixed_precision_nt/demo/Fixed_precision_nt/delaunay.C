#include <CGAL/basic.h>
#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <strstream.h>

#include <CGAL/Fixed_precision_nt.h>

#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h>   // to avoid a g++ problem
#include <CGAL/Point_2.h>

#include <CGAL/Triangulation_euclidean_traits_2.h>

#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/IO/Window_stream.h>

using namespace CGAL;


typedef Fixed_precision_nt coord_type;
static bool Fixed_precision_nt_init_result 
                =  Fixed_precision_nt::init(2000.0);

typedef Cartesian<coord_type>  Repclass;

typedef Point_2<Repclass>  Point_;

typedef Triangulation_euclidean_traits_2<Repclass> Traits_;
typedef Triangulation_vertex_base_2<Traits_> Vb;
typedef Triangulation_face_base_2<Traits_>  Fb;
typedef Triangulation_default_data_structure_2<Traits_,Vb,Fb> Tds;
typedef Delaunay_triangulation_2<Traits_,Tds>  Delaunay_;




int main(int argc, char* argv[])
{
    Delaunay_ D;
    Window_stream W(200,200); // physical window size

    W.init(-1,1,-1);   // logical window size
    W << BLUE;
    W.set_mode(leda_src_mode);
    W.set_node_width(3);
    W.display();
    
    cout <<endl<<endl<< "DELAUNAY TRIANGULATION"<<endl;
    
    while(1) {
      double x, y;
      int b = W.get_mouse(x,y);
      if ( ! ( (b == MOUSE_BUTTON(1)) ||
	       (b == MOUSE_BUTTON(2)) ||
	       (b == MOUSE_BUTTON(3)))) continue;
      if (b!=MOUSE_BUTTON(1)) break;
      D.insert( Point_(coord_type(x), coord_type(y)));
      W.clear();
      W << D;
    }
    return 0;
}


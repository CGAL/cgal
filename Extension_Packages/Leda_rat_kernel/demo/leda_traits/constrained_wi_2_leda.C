#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/IO/Window_stream.h>
#include <LEDA/window.h>
#include <LEDA/rat_window.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>


#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif


typedef CGAL::leda_rat_kernel_traits       K;

typedef K::Point_2     Point;
typedef K::Segment_2   Segment;

typedef CGAL::Triangulation_vertex_base_2<K>                   Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>         Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>            TDS;
typedef CGAL::Exact_intersections_tag                          Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS,Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>            CDTplus;

typedef CGAL::Window_stream  Window_stream;


void draw_constraints(Window_stream&  W, const CDTplus& tr)
{
  CDTplus::Subconstraint_iterator scit=tr.subconstraints_begin();
  for ( ; scit != tr.subconstraints_end(); scit++) {
    	W << CGAL::RED << Segment( scit->first.first->point(),
				   scit->first.second->point());
  }
}

void window_input(std::list<Segment> & list_contraintes,  
	     CDTplus & t,
	     Window_stream &W)
{
  std::cerr << "Enter segments  with the left button" << std::endl;
  std::cerr << "double clic for a single point " << std::endl;
  std::cerr << "Right button terminates input" << std::endl;
  

  Point p,q;
  Segment s;
  double x,y;

  while(1) {
    int b = W.read_mouse(x,y);
    if(b== MOUSE_BUTTON(3)) break;
    if(b == MOUSE_BUTTON(1)) {
      p = Point(leda_point(x,y));
      W << p;
      W >> q;
      s = Segment(p,q);

      W << s;
      list_contraintes.push_back(s);

      W.set_mode(leda_src_mode);

      t.insert(p,q);
      assert(t.is_valid());
    }
    W.clear();
    W << CGAL::BLUE << t;
    draw_constraints(W,t);
  }
}

int main()
{
   Window_stream W(400, 400); // physical window size
   W.init(-1, 1, 1);   // logical window size
   W << CGAL::BLUE;
   CGAL::cgalize( W);
   W.display();
    
   CDTplus t;
   std::list<Segment> list_contraintes;
   window_input(list_contraintes, t, W); 
   return 0;
}


#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) && !defined(CGAL_USE_CGAL_WINDOW)
int main(int argc, char* argv[])
{

  std::cout << "Sorry, this demo needs LEDA for visualisation. or CGAL_WINDOW";
  std::cout << std::endl;

  return 0;
}

#else
#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_drawing_mode CGAL::drawing_mode
#define leda_xor_mode     CGAL::xor_mode
#define leda_src_mode     CGAL::src_mode
#define leda_red          CGAL::red
#define leda_yellow       CGAL::yellow
#define leda_invisible    CGAL::invisible
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/intersections.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/IO/Window_stream.h>

typedef CGAL::Simple_cartesian<double>  K1;
typedef CGAL::Filtered_kernel<K1>       K2;
struct K : public K2 {};

typedef K::Point_2  Point;
typedef K::Segment_2 Segment;
typedef CGAL::Triangulation_vertex_base_2<K>              Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>       TDS;
typedef CGAL::Exact_predicates_tag                        Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS,Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT> CDTplus;

typedef CGAL::Window_stream  Window_stream;

void
any_button(CGAL::Window_stream &W)
{
    double x, y;
    std::cerr << "Press any button to continue" << std::endl;
    W.read_mouse(x,y);
}

void
draw_constraints(Window_stream&  W  , const CDTplus& tr)
{
  CDTplus::Subconstraint_iterator scit=tr.subconstraints_begin();
  for ( ; scit != tr.subconstraints_end(); scit++) {
    	W << CGAL::RED << Segment( scit->first.first->point(),
				   scit->first.second->point());
  }
}

void
window_input(std::list<Segment> & list_contraintes,  
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
      p = Point(x,y);
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

#endif // CGAL_USE_LEDA || CGAL_USE_CGAL_WINDOW

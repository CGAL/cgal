// examples/Pm_with_intersections/example2
// ---------------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>

// We use here double instead of leda_rational to enable compilation 
// without LEDA.
// This is not recommended generally.
// Read more in the README file or in the manual.
typedef double                                              NT;
typedef CGAL::Cartesian<NT>                                 R;
typedef CGAL::Arr_segment_exact_traits<R>                   Traits;

typedef Traits::Point                                       Point;
typedef Traits::X_curve                                     X_curve;

typedef CGAL::Pm_default_dcel<Traits>                       Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                     Planar_map_2;
typedef CGAL::Planar_map_with_intersections_2<Planar_map_2> Pmwx;

typedef Pmwx::Pmwx_change_notification Pmwx_change_notification;

using namespace std;

class My_notification : public Pmwx_change_notification 
{
public:
	
  My_notification()
  {i = 0;}
	  
  void add_edge(const  Traits::X_curve& cv, Planar_map::Halfedge_handle e, 
		bool left_to_right, bool overlap=false)
  {
    cout << "add_edge" << endl;
    i++;
  }
	
  void split_edge(Planar_map::Halfedge_handle orig_edge, 
		  Planar_map::Halfedge_handle new_edge,
		  const Traits::X_curve& c1, const Traits::X_curve& c2)
  {
    cout << "split_edge" << endl;
    i++;
  }
	
  void split_face(Planar_map::Face_handle orig_face, 
		  Planar_map::Face_handle new_face)
  {
    cout << "split_face" << endl;
  }

  void add_hole(Planar_map::Face_handle in_face, 
		Planar_map::Halfedge_handle new_hole)
  {
    cout << "add_hole" << endl;
  }
	
  int i;
};

int main() {
  
  Pmwx pm;
  My_notification notif;

  //insertion of the curves
  X_curve c1(Point(0,1),Point(1,0));
  X_curve c2(Point(0,0),Point(1,1));
  X_curve c3(Point(0,1),Point(1,1));

  cout << "inserting " << c1 << endl;
  pm.insert(c1, &notif);
  cout << "inserting " << c2 << endl;
  pm.insert(c2, &notif);
  cout << "inserting " << c3 << endl;
  pm.insert(c3, &notif);

  cout << "Total number of edges " << notif.i << endl;

  return 0;
}

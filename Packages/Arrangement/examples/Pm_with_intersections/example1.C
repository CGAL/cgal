// examples/Pm_with_intersections/example1
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
typedef double                                        NT;
typedef CGAL::Cartesian<NT>                           R;
typedef CGAL::Arr_segment_exact_traits<R>             Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;

typedef CGAL::Pm_default_dcel<Traits>                          Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                        Planar_map_2;
typedef CGAL::Planar_map_with_intersections_2<Planar_map_2>    Pmwx;

using namespace std;

int main() {
  
  Pmwx pm;
  X_curve cv1(Point(0,0),Point(1,1));
  X_curve cv2(Point(0,1),Point(1,0)); 

  //insertion of the curves
  cout << "Inserting the segments:" << endl;
  cout << cv1 << endl;
  pm.insert(cv1);
  cout << cv2 << endl << endl;
  pm.insert(cv2);
  
  //traversal of the curves
  cout << "Edges of the planar map:" << endl;

  Pmwx::Halfedge_iterator eit;
  for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) {
    cout << eit->source()->point();
    cout << " --- ";
    cout << eit->target()->point();
    cout << endl;
  }

  return 0;
}

//example7.C
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Arrangement_2.h>

// We use here double instead of leda_real to enable compilation without LEDA.
// This is not recommended generally.
// Read more in the README file or in the manual.
typedef double                                        NT;
typedef CGAL::Cartesian<NT>                            R;
typedef CGAL::Arr_segment_exact_traits<R>              Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

int global_int=0; // a global variable to keep track of when we inserted the edges.

struct Info { 
  Info() : d(global_int) {} 
  //initializes the private instance with the global variable
  
  int num() {return d;}
  int d;
};

struct My_base_node : public CGAL::Arr_base_node<Curve> {
  typedef CGAL::Arr_base_node<Curve> Base;
  My_base_node() : Base(), info() {}

  Info info; //will be initialized with the default ctr.
};


struct My_halfedge : public CGAL::Arr_2_halfedge_base<My_base_node> {
  typedef CGAL::Arr_2_halfedge_base<My_base_node> Base;
  My_halfedge() : Base(), info() {}

  Info info; //will be initialized with the default ctr.
};

typedef CGAL::Pm_dcel<CGAL::Arr_2_vertex_base< Point >,  
                      My_halfedge, 
                      CGAL::Arr_2_face_base >             Dcel;

typedef CGAL::Arrangement_2<Dcel,Traits,My_base_node >    Arr_2;

using namespace std;

int main() {
  Arr_2 arr;
  
  arr.set_update(false); //the insertions won't update the planar map.

  //insertion of the curves
  Arr_2::Curve_iterator cit=arr.insert(Curve(Point(0,0),Point(1,1)));
  cit=arr.insert(Curve(Point(0,1),Point(1,0))); 
  CGAL_assertion(arr.number_of_halfedges()==0);
 
  global_int=1; //change the global variable so we will see the changes
  //all Info structs will have 1 in them from now on.

  arr.set_update(true);
  CGAL_assertion(arr.number_of_halfedges()==8);

  //traversal of the curves
  Arr_2::Edge_iterator eit;
  for (cit=arr.curve_node_begin(); cit!=arr.curve_node_end(); ++cit) {
    cout << "\nCurve level info:\n" << cit->info.num() << endl ;
    cout << "Edge level info:\n";
    for (eit=cit->edges_begin(); eit!=cit->edges_end(); ++eit) {
      cout << eit->info.num() << endl ;
    }
  }
  return 0;
}

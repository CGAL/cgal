//example4.C
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

// We use here double instead of leda_real to enable compilation without LEDA.
// This is not recommended generally.
// Read more in the README file or in the manual.
typedef double                                        NT;
typedef CGAL::Cartesian<NT>                            R;
typedef CGAL::Arr_segment_exact_traits<R>              Traits;
typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;
typedef CGAL::Arr_base_node<Curve>                     Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>               Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >    Arr_2;

using namespace std;

//a simple functions that splits a segment into 2
void my_split_f(const Curve& cv, list<Curve>& l) {
  Point s=cv.source(); //uses the knowledge of the curve functions
  Point t=cv.target();
  Point m1=s+(t-s)/2.0;
  l.push_back(Curve(s,m1));
  l.push_back(Curve(m1,t));
}

typedef void (*SPLIT_FUNC)(const Curve& cv, list<Curve>& l);

int main() {
   vector<SPLIT_FUNC> func_vec;
   func_vec.push_back(&my_split_f);
   Arr_2 arr;

   //insertion with user-defined function
   Arr_2::Curve_iterator cit=arr.insert(Curve(Point(0,0),Point(6,6)),
                                     func_vec.begin(),func_vec.end());

   //regular insertion                                  
   cit=arr.insert(Curve(Point(0,4),Point(6,4))); 
  
   //traversal of the curves
   Arr_2::Edge_iterator eit;
   for (cit=arr.curve_node_begin(); cit!=arr.curve_node_end(); ++cit) {
      cout << "\nCurve level:\n" << cit->curve() << endl ;
      cout << "Edge level:\n";
      for (eit=cit->edges_begin(); eit!=cit->edges_end(); ++eit) {
         cout << eit->curve() << endl ;
      }
   }

   return 0;
}

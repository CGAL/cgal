//example2.C
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Arrangement_2.h>

// We use here double instead of leda_real to enable compilation without LEDA.
// This is not recommended generally.
// Read more in the README file or in the manual.
typedef CGAL::Cartesian<double>                            R;
typedef CGAL::Arr_segment_exact_traits<R>              Traits;

typedef Traits::Point                                 Point;
typedef Traits::Curve                                 Curve;

typedef CGAL::Arr_base_node<Curve>                     Base_node;
typedef CGAL::Arr_2_default_dcel<Traits> Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >    Arr_2;

using namespace std;

int main() {
   Arr_2 arr;

   //insertion of the curves
   Arr_2::Curve_iterator cit=arr.insert(Curve(Point(0,0),Point(2,2)));
   cit=arr.insert(Curve(Point(1,1),Point(3,3))); 
  
   //traversal of the halfedges
   Arr_2::Halfedge_const_iterator hit;
   for (hit=arr.halfedges_begin(); hit!=arr.halfedges_end(); ++hit,++hit) {
     //we skip the adjacent twin halfedge
     Arr_2::Overlap_const_circulator occ=hit->overlap_edges();
     int count=0;
     do {
       count++;
     } while (++occ!=hit->overlap_edges());

   if (count == 1) 
     cout << "Edge " << occ->curve() << " is covered by a single edge.\n";
   else
     cout << "Edge " << occ->curve() << " is covered by " << count << " edges.\n";

   }

   return 0;
}

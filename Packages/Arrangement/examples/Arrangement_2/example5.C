//examples/Arrangement_2/example5.C
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
typedef CGAL::Cartesian<NT>                           R;
typedef CGAL::Arr_segment_exact_traits<R>             Traits;
typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;
typedef CGAL::Arr_base_node<Curve>                    Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>              Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >   Arr_2;

using namespace std;

//a base class for split functions
struct Split_base {
  virtual void operator()(const Curve& cv, list<Curve>& l)=0;
};

struct Split_func : public Split_base {
  Split_func(double ratio) : r(ratio) {}
  void operator()(const Curve& cv, list<Curve>& l) {
     Point s=cv.source(); //uses the knowledge of the curve functions
     Point t=cv.target();
     Point m1=s+(t-s)/r;
     l.push_back(Curve(s,m1));
     l.push_back(Curve(m1,t));
   }     
private:
  double r;
};

int main() {
   std::vector<Split_base*> func_vec;
   func_vec.push_back(new Split_func(2.0));
   func_vec.push_back(new Split_func(3.0));
   Arr_2 arr;

   //insertion with user-defined function
   Arr_2::Curve_iterator cit=arr.insert(Curve(Point(0,0),Point(6,6)),
                                     func_vec.begin(),func_vec.end());

   CGAL_assertion(arr.number_of_halfedges()==8);                               

   CGAL_assertion(cit->number_of_sc_levels()==2);
   return 0;
}

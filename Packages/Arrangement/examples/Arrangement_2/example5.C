// examples/Arrangement_2/example5.C
// ---------------------------------

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_exact_traits<Kernel>          Traits;
typedef Traits::Point_2                                 Point;
typedef Traits::Curve_2                                 Curve;
typedef CGAL::Arr_base_node<Curve>                      Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node>      Arr_2;

// A base class for split functors
struct Split_base 
{
  virtual void operator()(const Curve& cv, std::list<Curve>& l)=0;
};

// A user-defined insertion functor
struct Split_func : public Split_base 
{
  Split_func(double ratio) : r(ratio) {}
  void operator()(const Curve & cv, std::list<Curve> & l) 
  {
     Point s=cv.source(); // Uses the knowledge of the curve functions
     Point t=cv.target();
     Point m1 = s + (t - s) / r;
     l.push_back(Curve(s, m1));
     l.push_back(Curve(m1, t));

   }     

  virtual ~Split_func(){};
private:
  NT r;
};

int main() 
{
  // Prepare a vector of pointers to the functor base class 
  std::vector<Split_base*> func_vec;

  // Create 2 functors
  Split_func Sf1(2.0), Sf2(3.0);

  func_vec.push_back(&Sf1);
  func_vec.push_back(&Sf2);

  Arr_2 arr;

  // Insertion with user-defined functor
  Arr_2::Curve_iterator cit=arr.insert(Curve(Point(0, 0),Point(6, 6)),
                                       func_vec.begin(), func_vec.end());

  CGAL_assertion(arr.number_of_halfedges() == 8);
  CGAL_assertion(cit->number_of_sc_levels() == 2);

  return 0;
}

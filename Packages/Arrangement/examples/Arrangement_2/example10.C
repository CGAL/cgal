//example10.C

// Define shorter names to please linker (g++/egcs)
#define Arrangement_2 Ar
#define _In_place_list_iterator IPLI
#define Cartesian CRTS
#define Arr_polyline_traits APT

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_polyline_traits.h>

//#include <CGAL/leda_real.h>
//typedef leda_real                                     NT;

// We use here double instead of leda_real to enable compilation without LEDA.
// This is not recommended generally.
// Read more in the README file or in the manual.
typedef double                                        NT;
typedef CGAL::Cartesian<NT>                           Rep;
typedef CGAL::Arr_polyline_traits<Rep>                Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Arr_base_node<Curve>                    Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>              Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >   Arr_2;

using namespace std;

int main() {
   Arr_2 arr;
   Curve in_curve;

   // curve #1, not x monotone
   in_curve.push_back(Point(  0,  0));
   in_curve.push_back(Point( 10, 10));
   in_curve.push_back(Point(  0, 20));
   arr.insert(in_curve); 

   // curve #2, x monotone
   in_curve.clear();
   in_curve.push_back(Point(100,  0));
   in_curve.push_back(Point(150, 50));
   in_curve.push_back(Point(200,  0));
   arr.insert(in_curve);

   // curve #1 is broken into two edges. Point (10,10) turns into a vertex.
   Arr_2::Locate_type lt;
   arr.locate(Point(10, 10), lt);
   CGAL_assertion(lt == Arr_2::VERTEX);

   return 0;
}

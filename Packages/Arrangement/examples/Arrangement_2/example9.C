//examples/Arrangement_2/example9.C

// Define shorter names to please linker (g++/egcs)
#define Arrangement_2 Ar
#define _In_place_list_iterator IPLI
#define Cartesian CRTS
#define Arr_polyline_traits APT

#include <CGAL/Cartesian.h>

#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_polyline_traits.h>
#include <CGAL/Arrangement_2.h>

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

   // we insert two intersecting squares
   in_curve.push_back(Point(  0, 50));
   in_curve.push_back(Point( 50, 50));
   in_curve.push_back(Point( 50,  0));
   in_curve.push_back(Point(  0,  0));
   in_curve.push_back(Point(  0, 50));
   arr.insert(in_curve);

   in_curve.clear();
   in_curve.push_back(Point( 25, 75));
   in_curve.push_back(Point( 75, 75));
   in_curve.push_back(Point( 75, 25));
   in_curve.push_back(Point( 25, 25));
   in_curve.push_back(Point( 25, 75));
   arr.insert(in_curve); 

   // upward vertical ray shooting
   // the edge <25,50>-<50,50> is supposed to be have been created
   Arr_2::Locate_type lt;
   Arr_2::Halfedge_handle e=arr.vertical_ray_shoot(Point(30,30),lt,true);

   CGAL_assertion(e->source()->point()==Point(50,50)); 
   CGAL_assertion(e->target()->point()==Point(25,50));

   // we expect <50,25> to be an intersection (of the polylines)
   e = arr.locate(Point(50,25), lt);
   CGAL_assertion(lt == Arr_2::VERTEX);

   return 0;
}

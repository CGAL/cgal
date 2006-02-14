// provide 3d kernel traits ...
#define CGAL_NO_POSTCONDITIONS  
#define CGAL_NO_PRECONDITIONS  
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3

#include <CGAL/Cartesian.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>

// include spezializations for LEDA d3_rat_point's
// of the projective traits ...
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/Convex_hull_projective_xy_traits_leda_rat_2.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/Convex_hull_projective_xz_traits_leda_rat_2.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/Convex_hull_projective_yz_traits_leda_rat_2.h>

#include <CGAL/Convex_hull_traits_3.h>

namespace CGAL {
template <>
class Max_coordinate_3<leda_rat_vector> 
{
public:

    int operator()(const leda_rat_vector& v)
    {
      if (CGAL_NTS abs(v.xcoord()) >= CGAL_NTS abs(v.ycoord()))
      {
         if (CGAL_NTS abs(v.xcoord()) >= CGAL_NTS abs(v.zcoord())) return 0;
         return 2;
      }
      else
      {
         if (CGAL_NTS abs(v.ycoord()) >= CGAL_NTS abs(v.zcoord())) return 1;
         return 2;
      }
    }
};
}

#include <CGAL/convex_hull_3.h>

#include <iostream>
#include <list>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

using std::cin;
using std::cout;

typedef CGAL::leda_rat_kernel_traits              K;

// Probleme : der instanziierte Convex_hull_projective_xy_traits_2.h
// braucht x() und y() Memberfunktionen auf den Punkten ... 

typedef CGAL::Convex_hull_traits_3<K>             Traits;
typedef Traits::Polyhedron_3                      Polyhedron_3;
typedef K::Segment_3                              Segment_3;
typedef K::Point_3                                Point_3;

// generator ...
void generate_input(int n, int input, int maxc, leda_list<leda_d3_rat_point>& L)
{ 
  L.clear();

  switch (input) {
   case 0: random_points_in_cube(n,3*maxc/4,L); break;
   case 1: random_points_in_ball(n,maxc,L); break;
   case 2: random_points_in_square(n,maxc,L); break;
   case 3: random_points_on_paraboloid(n,maxc,L); break;
   case 4: lattice_points(n,3*maxc/4,L); break;
   case 5: random_d3_rat_points_on_sphere(n,maxc,L);  break;
   case 6: random_points_on_segment(n,maxc,L); break;
  }
}

int main()
{
  leda_list<Point_3> L2;
  
  int numb;
  int gen;
 
  cout << "#Points   :"; cin >> numb;
  cout << "0 - in cube; 1 - in ball; 2 - in square; 3 - on paraboloid;\n";
  cout << "4 - lattice points; 5 - on sphere; 6 - on segment\n";
  cout << "Type of point generation:"; 
  cin >> gen;
 
  generate_input(numb, gen, 1000, L2);

  
  // define object to hold convex hull 
  CGAL::Object ch_object;
  K            kernel;
  Traits       tr;

  // compute convex hull 
  CGAL::convex_hull_3(L2.begin(), L2.end(), ch_object, tr);

  // determine what kind of object it is
  Segment_3 segment;
  Polyhedron_3 polyhedron;
  if ( CGAL::assign(segment, ch_object) )
     std::cout << "convex hull is a segment " << std::endl;
  else if ( CGAL::assign (polyhedron, ch_object) )
     std::cout << "convex hull is a polyhedron " << std::endl;
  else
     std::cout << "convex hull error!" << std::endl;

  return 0;
}

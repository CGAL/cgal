// provide 3d kernel traits ...
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3

#include <CGAL/basic.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

#include <LEDA/integer.h>
#include <LEDA/misc.h>

#include <iostream>
#include <list>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

using std::cin;
using std::cout;

typedef CGAL::leda_rat_kernel_traits                     K;

typedef CGAL::Triangulation_vertex_base_3<K>             Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
typedef CGAL::Triangulation_cell_base_3<void>            Cb;
typedef CGAL::Triangulation_data_structure_3<Vbh,Cb>     Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds>            DT;
typedef CGAL::Triangulation_hierarchy_3<DT>              DH;

typedef K::Point_3                                       Point;

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
  DT T;
  DH HT; // hierarchy
  // insertion of points
  int x,y,z;
  leda_list<Point> L2;
  
  int numb;
  int gen;
 
  cout << "#Points   :"; cin >> numb;
  cout << "0 - in cube; 1 - in ball; 2 - in square; 3 - on paraboloid;\n";
  cout << "4 - lattice points; 5 - on sphere; 6 - on segment\n";
  cout << "Type of point generation:"; 
  cin >> gen;
 
  generate_input(numb, gen, 1000, L2);
 
  cout << "We compute the Delaunay triangulation ...\n";
  float tm = used_time(); 	  
  T.insert(L2.begin(), L2.end());
  cout << "time: " << used_time(tm) << " seconds.\n";
  
  // check the result ...
  bool valid = T.is_valid(true);
  
  if (valid) std::cout << "valid result !\n";
  else std::cout << "NON - valid result !\n";
 
  cout << "We compute the Delaunay triangulation using the hierarchy...\n";
  tm = used_time(); 	  
  HT.insert(L2.begin(), L2.end());
  cout << "time: " << used_time(tm) << " seconds.\n";
  
  // check the result ...
  valid = HT.is_valid(true);
  
  if (valid) std::cout << "valid result !\n";
  else std::cout << "NON - valid result !\n"; 
  
  return 0;
}

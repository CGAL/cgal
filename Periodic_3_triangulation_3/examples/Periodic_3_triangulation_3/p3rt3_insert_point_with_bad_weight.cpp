#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>

#include <iostream>


typedef CGAL::Epick K;
typedef K::FT FT;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   RT;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<RT> Traits;

/* If remove() isn't called in our program, we can use a triangulation data structure more appropriate
 * which saves some memory resources.
 */
typedef CGAL::Triangulation_vertex_base_3<Traits, CGAL::Periodic_3_triangulation_ds_vertex_base_3<> > Vb;
typedef CGAL::Triangulation_cell_base_3<Traits, CGAL::Periodic_3_triangulation_ds_cell_base_3<> > Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;

typedef CGAL::Periodic_3_regular_triangulation_3<Traits, Tds>    P3RT3;

typedef Traits::Iso_cuboid_3 Iso_cuboid;
typedef Traits::Weighted_point Weighted_point;
typedef Traits::Bare_point Bare_point;


int main ()
{
  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

//  p3rt3.insert(Weighted_point(Bare_point(0.5,0.5,0.5),1.)); // Here, we insert a point with a bad weight.

//   In debug mode, if we uncomment the previous instruction, the program displays the following error message :
//   terminate called after throwing an instance of 'CGAL::Precondition_exception'
//     what():  CGAL ERROR: precondition violation!
//   Expr: point.weight() < ( FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin()) )
//   [...]
//   Explanation: point.weight() < 1/64 * domain_size * domain_size

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}

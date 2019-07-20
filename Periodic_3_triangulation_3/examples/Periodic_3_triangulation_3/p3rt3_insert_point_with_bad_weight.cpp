#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>

#include <iostream>

typedef CGAL::Epick                                         K;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>  Gt;

/* If remove() isn't called in our program, we can use a triangulation data structure more appropriate
 * which saves some memory resources.
 */
typedef CGAL::Regular_triangulation_vertex_base_3<Gt,
          CGAL::Periodic_3_triangulation_ds_vertex_base_3<> >    Vb;
typedef CGAL::Regular_triangulation_cell_base_3<Gt,
          CGAL::Periodic_3_triangulation_ds_cell_base_3<> >      Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>             Tds;

typedef CGAL::Periodic_3_regular_triangulation_3<Gt, Tds>        P3RT3;

typedef Gt::Iso_cuboid_3                                         Iso_cuboid;
typedef Gt::Weighted_point_3                                     Weighted_point_3;
typedef Gt::Point_3                                              Point_3;

int main ()
{
  P3RT3 p3rt3(P3RT3::Iso_cuboid(0,0,0, 1,1,1));

  // Here, we insert a point with a bad weight.
//  p3rt3.insert(Weighted_point_3(Point_3(0.5,0.5,0.5), 1.));

//   In debug mode, if we uncomment the previous instruction, the program displays the following error message :
//   terminate called after throwing an instance of 'CGAL::Precondition_exception'
//     what():  CGAL ERROR: precondition violation!
//   Expr: point.weight() < ( FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin()) )
//   [...]
//   Explanation: point.weight() < 1/64 * domain_size * domain_size

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}

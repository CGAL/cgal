#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT                                               Weight;
typedef K::Point_3                                          Point;
typedef K::Weighted_point_3                                 Weighted_point;

typedef CGAL::Regular_triangulation_vertex_base_3<K>        Vb0;

typedef CGAL::Triangulation_vertex_base_with_info_3<int, K, Vb0> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<K>          Cb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb>         Tds;

typedef CGAL::Regular_triangulation_3<K, Tds>               Rt;

int main()
{
  Weighted_point wp(CGAL::ORIGIN,1);

  Rt rt;

  rt.insert(wp);

  assert( rt.is_valid() );

  return 0;
}

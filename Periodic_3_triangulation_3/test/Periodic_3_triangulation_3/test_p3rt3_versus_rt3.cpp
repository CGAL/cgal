#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Regular_triangulation_3.h>

#include <CGAL/periodic_3_triangulation_3_io.h>
#include <CGAL/IO/Triangulation_off_ostream_3.h>
#include <CGAL/Timer.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

// A basic test to check that P3RT3 and RT3 produces the same regular triangulations
// In the case of RT3, a fake periodicity is obtained by adding 26 copies
// of the same input point set (similarly to how we copy in P3RT3 when it cannot
// yet be converted to 1 sheet) around the cube.

typedef CGAL::Exact_predicates_inexact_constructions_kernel       Epick;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<Epick>    PRTT_Inexact;

typedef CGAL::Periodic_3_regular_triangulation_3<PRTT_Inexact>    P3RT3;

typedef P3RT3::Iso_cuboid                                         Iso_cuboid;
typedef P3RT3::Offset                                             Offset;
typedef P3RT3::Weighted_point                                     Weighted_point;

// The info() is a unique color per point + a boolean to mark those in the "fake"
// instances (the cuboid surrounded by 26 copies)
typedef CGAL::Regular_triangulation_vertex_base_3<Epick>          Vb0;
typedef CGAL::Triangulation_vertex_base_with_info_3<
                std::pair<int, bool>, Epick, Vb0>                 Vb;
typedef CGAL::Regular_triangulation_cell_base_3<Epick>            Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>              Tds;

typedef CGAL::Regular_triangulation_3<Epick, Tds>                 RT3;
typedef RT3::Vertex_handle                                        Vertex_handle;
typedef RT3::Finite_vertices_iterator                             Finite_vertices_iterator;
typedef RT3::Finite_cells_iterator                                Finite_cells_iterator;

int main (int, char**)
{
  CGAL::Timer t;
  t.start();
  Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  P3RT3 p3rt3(iso_cuboid);

  std::ifstream input_stream("data/p3rt3_point_set__s7_n800");

  std::vector<Weighted_point> points_for_p3rt3;
  points_for_p3rt3.reserve(800);

  while(points_for_p3rt3.size() != 800)
  {
    Weighted_point p;
    input_stream >> p;
    points_for_p3rt3.push_back(p);
  }
  p3rt3.insert(points_for_p3rt3.begin(), points_for_p3rt3.end(), true /*large point set*/);
  std::cout << "--- built periodic regular triangulation" << std::endl;

  // Non-periodic
  RT3 rt3;

  int id = 0;
  for(int l=0; l<800; ++l)
  {
    for(int i=-1; i<2; i++)
    {
      for(int j=-1; j<2; j++)
      {
        for(int k=-1; k<2; k++)
        {
          Vertex_handle v = rt3.insert(
                              p3rt3.geom_traits().construct_weighted_point_3_object()(
                                points_for_p3rt3[l], Offset(i, j, k)));

          if(v != Vertex_handle())
            v->info() = std::make_pair(id /*unique id*/,
                                       (i==0 && j==0 && k==0) /*canonical instance?*/);
        }
      }
    }
    id++;
  }
  std::cout << "--- built regular triangulation" << std::endl;

//  std::ofstream out_p3rt3("out_p3rt3.off");
//  std::ofstream out_rt3("out_rt3.off");
//  CGAL::write_triangulation_to_off(out_p3rt3, p3rt3);
//  CGAL::export_triangulation_3_to_off(out_rt3, rt3);

  // ---------------------------------------------------------------------------

  // compare vertices
  std::set<int> unique_vertices_from_rt3;

  Finite_vertices_iterator vit = rt3.finite_vertices_begin();
  for(; vit!=rt3.finite_vertices_end(); ++vit)
  {
    if(vit->info().second) // is in the canonical instance
      unique_vertices_from_rt3.insert(vit->info().first);
  }

  std::cout << unique_vertices_from_rt3.size() << " unique vertices (rt3)" << std::endl;
  std::cout << p3rt3.number_of_vertices() << " unique vertices (p3rt3)" << std::endl;
  assert(unique_vertices_from_rt3.size() == p3rt3.number_of_vertices());

  // compare cells
  std::set<std::set<int> > unique_cells_from_rt3;

  Finite_cells_iterator cit = rt3.finite_cells_begin();
  for(; cit!=rt3.finite_cells_end(); ++cit)
  {
    std::set<int> cell;
    bool has_a_point_in_canonical_domain = false;

    for(int i=0; i<4; i++)
    {
      Vertex_handle v = cit->vertex(i);
      cell.insert(v->info().first);

      if(v->info().second)
        has_a_point_in_canonical_domain = true;
    }

    if(has_a_point_in_canonical_domain)
      unique_cells_from_rt3.insert(cell);
  }

  std::cout << unique_cells_from_rt3.size() << " unique cells (rt3)" << std::endl;
  std::cout << p3rt3.number_of_cells() << " unique cells (p3rt3)" << std::endl;
  assert(unique_cells_from_rt3.size() == p3rt3.number_of_cells());

  std::cout << t.time() << " sec." << std::endl;
  std::cout << "EXIT SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#define CGAL_P2T2_USE_COMPACT_OFFSET

#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/periodic_2_triangulation_2_io.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/point_generators_2.h>

#include <CGAL/Timer.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

// A basic test to check that P2DT2 and DT2 produces the same Delaunay triangulations
// In the case of DT2, a fake periodifity is obtained by addding 26 copies
// of the same input point set (similarly to how we copy in P2DT2 when it cannot
// yet be converted to 1 sheet) around the cube.

typedef CGAL::Exact_predicates_inexact_constructions_kernel       Epick;
typedef Epick::FT                                                 FT;

typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<Epick>   Traits;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<Traits>         P2DT2;
typedef P2DT2::Domain                                             Domain;
typedef P2DT2::Offset                                             Offset;
typedef P2DT2::Point                                              Point;

// The info() is a unique color per point + a boolean to mark those in the "fake"
// instances (the cuboid surrounded by 26 copies)
typedef CGAL::Triangulation_vertex_base_2<Epick>                  Vb0;
typedef CGAL::Triangulation_vertex_base_with_info_2<
                std::pair<int, bool>, Epick, Vb0>                 Vb;
typedef CGAL::Triangulation_face_base_2<Epick>                    Cb;
typedef CGAL::Triangulation_data_structure_2<Vb, Cb>              Tds;
typedef CGAL::Delaunay_triangulation_2<Epick, Tds>                DT2;
typedef DT2::Vertex_handle                                        Vertex_handle;
typedef DT2::Finite_vertices_iterator                             Finite_vertices_iterator;
typedef DT2::Finite_faces_iterator                                Finite_faces_iterator;

typedef CGAL::Creator_uniform_2<FT, Point>                        Creator;
typedef CGAL::Random_points_in_square_2<Point, Creator>           Point_generator;

int main (int argc, char** argv)
{
  CGAL::Timer t;
  t.start();
  Domain domain(-0.5, -0.5, 0.5, 0.5);
  P2DT2 p2dt2(domain);

  std::size_t number_of_points = (argc > 1) ? std::atoi(argv[1]) : 800;
  std::vector<Point> points_for_p2dt2;
  points_for_p2dt2.reserve(number_of_points);

  int random_seed = (argc > 2) ? std::atoi(argv[2]) : CGAL::get_default_random().get_int(0, (1 << 30));
  std::cout << "random seed: " << random_seed << std::endl;

  CGAL::Random rnd(random_seed);
  std::copy_n(Point_generator(0.5, rnd), number_of_points, back_inserter(points_for_p2dt2));

  p2dt2.insert(points_for_p2dt2.begin(), points_for_p2dt2.end(), false /*large point set*/);

  std::cout << "--- built periodic Delaunay triangulation" << std::endl;
  std::cout << p2dt2.number_of_vertices() << " " << p2dt2.number_of_faces() << std::endl;
  std::cout << "1 cover?: " << p2dt2.is_1_cover() << std::endl;

  // Non-periodic
  DT2 dt2;

  int id = 0;
  for(std::size_t l=0; l<number_of_points; ++l)
  {
    for(int i=-1; i<2; ++i)
    {
      for(int j=-1; j<2; ++j)
      {
        Vertex_handle v = dt2.insert(
                            p2dt2.geom_traits().construct_point_2_object()(
                              points_for_p2dt2[l], Offset(i, j)));

        if(v != Vertex_handle())
          v->info() = std::make_pair(id /*unique id*/, (i==0 && j==0) /*canonical instance?*/);
      }
    }
    ++id;
  }
  std::cout << "--- built Delaunay triangulation" << std::endl;

  std::ofstream out_p2dt2("out_p2dt2.off");
  std::ofstream out_dt2("out_dt2.off");
  CGAL::write_P2T2_to_OFF(out_p2dt2, p2dt2);
  CGAL::draw_t2(out_dt2, dt2);

  // ---------------------------------------------------------------------------

  // compare vertices
  std::set<int> unique_vertices_from_dt2;

  Finite_vertices_iterator vit = dt2.finite_vertices_begin();
  for(; vit!=dt2.finite_vertices_end(); ++vit)
  {
    if(vit->info().second) // is in the canonical instance
      unique_vertices_from_dt2.insert(vit->info().first);
  }

  std::cout << unique_vertices_from_dt2.size() << " unique vertices (dt2)" << std::endl;
  std::cout << p2dt2.number_of_vertices() << " unique vertices (p2dt2)" << std::endl;
  CGAL_assertion(unique_vertices_from_dt2.size() == p2dt2.number_of_vertices());

  // compare faces
  std::set<std::map<Point, int> > unique_faces_from_dt2; // just to sort points lexicographically

  Finite_faces_iterator fit = dt2.finite_faces_begin();
  for(; fit!=dt2.finite_faces_end(); ++fit)
  {
    std::map<Point, int> face;
    for(int i=0; i<3; ++i)
    {
      Vertex_handle v = fit->vertex(i);
      face.emplace(v->point(), i);
    }

    Vertex_handle lowest_v = fit->vertex(face.begin()->second);
    if(lowest_v->info().second)
      unique_faces_from_dt2.insert(face);
  }

  std::cout << unique_faces_from_dt2.size() << " unique faces (dt2)" << std::endl;
  std::cout << p2dt2.number_of_faces() << " unique faces (p2dt2)" << std::endl;
  CGAL_assertion(unique_faces_from_dt2.size() == p2dt2.number_of_faces());

  std::cout << t.time() << " sec." << std::endl;
  std::cout << "EXIT SUCCESS" << std::endl;

  return EXIT_SUCCESS;
}

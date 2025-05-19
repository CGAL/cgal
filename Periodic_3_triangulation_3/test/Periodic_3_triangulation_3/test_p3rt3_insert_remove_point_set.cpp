//#define REALLY_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Timer.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <set>
#include <string>

typedef CGAL::Periodic_3_regular_triangulation_traits_3<CGAL::Epeck> Traits_Epeck;

typedef CGAL::Periodic_3_regular_triangulation_traits_3<CGAL::Epick> Traits_Epick;

template <class Kernel>
class Tests
{
public:
  typedef Kernel                                                K;
  typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>    Traits;

  typedef CGAL::Periodic_3_regular_triangulation_3<Traits>      P3RT3;
  typedef typename P3RT3::Vertex_handle                         Vertex_handle;
  typedef typename P3RT3::Facet                                 Facet;
  typedef typename P3RT3::Cell                                  Cell;
  typedef typename P3RT3::Cell_handle                           Cell_handle;
  typedef typename P3RT3::Vertex_iterator                       Vertex_iterator;
  typedef typename P3RT3::Cell_iterator                         Cell_iterator;
  typedef typename P3RT3::Segment                               Segment;
  typedef typename P3RT3::Triangle                              Triangle;
  typedef typename P3RT3::Tetrahedron                           Tetrahedron;
  typedef typename P3RT3::Periodic_tetrahedron                  Periodic_tetrahedron;
  typedef typename P3RT3::Locate_type                           Locate_type;
  typedef typename P3RT3::Offset                                Offset;

  typedef typename Traits::Weighted_point_3                     Weighted_point_3;
  typedef typename Traits::Point_3                              Point_3;
  typedef typename Traits::Iso_cuboid_3                         Iso_cuboid;

  static void test_insert_rnd_then_remove_all (unsigned pt_count)
  {
    typedef CGAL::Creator_uniform_3<double, Point_3>  Creator;
    CGAL::Random_points_in_cube_3<Point_3, Creator> in_cube(0.5);

    Iso_cuboid iso_cuboid(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
    P3RT3 p3rt3(iso_cuboid);

    std::vector<Weighted_point_3> insert_set;
    insert_set.reserve(pt_count);
    std::vector<Weighted_point_3> remove_set;
    remove_set.reserve(pt_count);

    std::cout << "-- insert" << std::endl;
    for (unsigned cnt = 1; cnt <= pt_count; ++cnt)
    {
      Weighted_point_3 p(*in_cube++, CGAL::get_default_random().get_double(0., 0.015625));

      std::size_t hidden_point_count = 0;
      for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
        hidden_point_count += iter->hidden_points().size();

      Vertex_handle vh = p3rt3.insert(p);

      std::size_t hidden_point_count_2 = 0;
      for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
        hidden_point_count_2 += iter->hidden_points().size();
      assert(hidden_point_count <= hidden_point_count_2);
      assert(hidden_point_count_2 + p3rt3.number_of_vertices() == cnt);
#ifdef REALLY_VERBOSE
      std::cout << cnt << " - p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;
#endif
      if (vh == Vertex_handle())
      {
        if (find(insert_set.begin(), insert_set.end(), p) == insert_set.end())
          insert_set.push_back(p);
      }
      else
        insert_set.push_back(p);
    }

    assert(p3rt3.is_valid());

    std::cout << "-- remove" << std::endl;
    unsigned cnt = 1;
    for (; p3rt3.number_of_vertices() != 0; ++cnt)
    {
      Vertex_iterator iter = p3rt3.vertices_begin();
      for (int j = CGAL::get_default_random().get_int(0, static_cast<int>(p3rt3.number_of_vertices())); j; --j)
        ++iter;
#ifdef REALLY_VERBOSE
      std::cout << cnt << " : " << iter->point() << std::endl;
#endif
      remove_set.push_back(iter->point());
      p3rt3.remove(iter);
#ifdef REALLY_VERBOSE
      std::cout << "    p3rt3.number_of_vertices() : " << p3rt3.number_of_vertices() << "  .number_of_stored_vertices : " << p3rt3.number_of_stored_vertices() << std::endl;
#endif
      std::size_t hidden_point_count = 0;
      for (Cell_iterator iter = p3rt3.cells_begin(), end_iter = p3rt3.cells_end(); iter != end_iter; ++iter)
        hidden_point_count += iter->hidden_points().size();

      assert(hidden_point_count + cnt + p3rt3.number_of_vertices() == insert_set.size());
    }

    std::sort(insert_set.begin(), insert_set.end());
    std::sort(remove_set.begin(), remove_set.end());
    assert(insert_set == remove_set);

    assert(p3rt3.is_valid());
    assert(p3rt3.number_of_stored_vertices() == 0);
  }

  static void test ()
  {
    test_insert_rnd_then_remove_all(800);
  }
};

int main (int, char**)
{
  CGAL::Timer t;
  t.start();
  std::cout << "TESTING ..." << std::endl;

  CGAL::Set_ieee_double_precision pfr;

  std::cout << "Epeck ..." << std::endl;
  Tests<CGAL::Epeck>::test();
  std::cout << "Epick ..." << std::endl;
  Tests<CGAL::Epick>::test();

  std::cout << t.time() << " sec." << std::endl;
  std::cout << "EXIT SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}

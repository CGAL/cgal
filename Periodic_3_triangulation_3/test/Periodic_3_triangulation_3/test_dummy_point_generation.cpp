#define CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>

#include <CGAL/IO/OFF.h>
#include <CGAL/Random.h>

#include <array>
#include <iostream>

// - Use a lattice rather than a grid to have a Delaunay separation because constructions are inexact
// - Use a reduced lattice as to use the 75 sufficient offsets giving the periodic star

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Iso_cuboid = K::Iso_cuboid_3;

using P3DTT3 = CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>;
using P3DT3 = CGAL::Periodic_3_Delaunay_triangulation_3<P3DTT3>;

using P3RTT3 = CGAL::Periodic_3_regular_triangulation_traits_3<K>;
using P3RT3 = CGAL::Periodic_3_regular_triangulation_3<P3RTT3>;

template <typename CellRangeIterator, typename Tr>
void dump_cells(const std::string filename,
                CellRangeIterator first,
                CellRangeIterator beyond,
                const Tr& tr,
                const Vector_3& translation)
{
  using Cell_handle = typename Tr::Cell_handle;

  auto cp = tr.geom_traits().construct_point_3_object();

  std::vector<Point_3> points;
  std::unordered_map<Point_3, std::size_t> points_to_ids;
  std::vector<std::vector<std::size_t> > faces;

  std::size_t pid = 0;
  auto add_point = [&](const Point_3& p) -> std::size_t
  {
    auto insertion_result = points_to_ids.emplace(p, pid);
    if(insertion_result.second)
    {
      points.push_back(p);
      ++pid;
    }

    return insertion_result.first->second;
  };

  for(auto cell_it = first; cell_it != beyond; ++cell_it)
  {
    const Cell_handle ch = *cell_it;

    const std::size_t v0_id = add_point(cp(tr.point(ch, 0)) + translation);
    const std::size_t v1_id = add_point(cp(tr.point(ch, 1)) + translation);
    const std::size_t v2_id = add_point(cp(tr.point(ch, 2)) + translation);
    const std::size_t v3_id = add_point(cp(tr.point(ch, 3)) + translation);

#if 1
    faces.push_back(std::vector<std::size_t>{v0_id, v1_id, v2_id});
    faces.push_back(std::vector<std::size_t>{v0_id, v2_id, v3_id});
    faces.push_back(std::vector<std::size_t>{v0_id, v3_id, v1_id});
    faces.push_back(std::vector<std::size_t>{v1_id, v2_id, v3_id});
#else
    faces.push_back(std::vector<std::size_t>{v0_id, v1_id, v2_id, v3_id});
#endif
  }

  CGAL::IO::write_OFF(filename, points, faces, CGAL::parameters::stream_precision(17));
}

template <typename Tr>
void dump(const std::string filename,
          const Tr& tr)
{
  typedef typename Tr::Cell_handle Cell_handle;
  std::vector<Cell_handle> chs;
  for(auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
    chs.push_back(cit);

  dump_cells(filename, std::cbegin(chs), std::cend(chs), tr, CGAL::NULL_VECTOR);
}

template <typename Tr>
void test(FT x_span, FT y_span, FT z_span)
{
  Tr p3t3;
  p3t3.set_domain(Iso_cuboid(0,0,0, x_span,y_span,z_span));
  p3t3.insert_generic_dummy_points();

  // Abuse the multi cover function to check if cells have a small-enough orthoradius
  p3t3.update_cover_data_after_converting_to_27_sheeted_covering();
  if(!p3t3.can_be_converted_to_1_sheet())
  {
    std::cerr << "Error: dummy points do not create a 1-cover" << std::endl;
    assert(false);
  }

//  dump("p3t3.off", p3t3);
}

int main (int argc, char** argv)
{
  CGAL::Random rnd = CGAL::get_default_random();
  std::cout << "seed: " << rnd.get_seed() << std::endl;

  if(argc != 1)
  {
    assert(argc == 4);
    test<P3DT3>(std::stod(argv[1]), std::stod(argv[2]), std::stod(argv[3]));
    std::cout << "OK" << std::endl;
    return EXIT_SUCCESS;
  }

  for(int i=0; i<10; ++i)
  {
    test<P3DT3>(rnd.get_double(0.01, 1.),
                rnd.get_double(0.01, 1.),
                rnd.get_double(0.01, 1.));
  }

  for(int i=1; i<5; ++i)
  {
    for(int j=1; j<5; ++j)
    {
      for(int k=1; k<5; ++k)
      {
#if (CGAL_P3T3_DUMMY_GENERATION_DEBUG_VERBOSITY >= 1)
        std::cout << "Test " << i << " " << j << " " << k << std::endl;
#endif
        test<P3DT3>(i,j,k);
        test<P3RT3>(i,j,k);
      }
    }
  }

  std::cout << "OK" << std::endl;
  return EXIT_SUCCESS;
}

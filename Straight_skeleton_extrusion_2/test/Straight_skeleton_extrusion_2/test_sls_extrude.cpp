#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/extrude_skeleton.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/IO/polygon_mesh_io.h>

namespace PMP = ::CGAL::Polygon_mesh_processing;

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using EPECK_w_SQRT = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;

template <typename FT>
bool are_equal(const FT t1, const FT t2,
               const FT delta = std::numeric_limits<FT>::epsilon(),
               const bool verbose = true)
{
  if(verbose)
    std::cerr << "Comparing " << t1 << " and " << t2 << " with authorized delta " << delta << std::endl;

  const FT diff = CGAL::abs(t1 - t2);
  if(verbose)
    std::cerr << "Diff: " << CGAL::abs(t1 - t2) << " vs eps: " <<  delta * (CGAL::abs(t1) + CGAL::abs(t2)) << std::endl;

  if(diff > delta && diff > delta * (CGAL::abs(t1) + CGAL::abs(t2)))
  {
    if(verbose)
    {
      std::cerr << "Approximate comparison failed (t1|t2): got " << t1 << " but expected " << t2 << std::endl;
    }
    return false;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename K>
bool read_dat_polygon(const char* filename,
                      CGAL::Polygon_with_holes_2<K>& p)
{
  using Point_2 = typename K::Point_2;
  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Error: Could not read " << filename << std::endl;
    return false;
  }

  bool is_number_of_CC_in_input = false;
  if(CGAL::IO::internal::get_file_extension(filename) == "poly")
  {
    is_number_of_CC_in_input = true;
  }

  std::vector<Polygon_2> polys;

  auto read_polygon = [&in, &polys](int i) -> void
  {
    std::vector<Point_2> poly;

    int v_count = 0;
    in >> v_count;
    for(int j=0; j<v_count && in; ++j)
    {
      double x = 0., y = 0.;
      in >> x >> y;
      poly.push_back(Point_2(x, y));
    }

    if(poly.size() >= 3)
    {
      bool is_simple = CGAL::is_simple_2(poly.begin(), poly.end(), K());
      if(!is_simple)
        std::cerr << "Warning: input polygon not simple (hopefully it is strictly simple...)" << std::endl;

      CGAL::Orientation expected = (i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE);

      const double area = CGAL::to_double(CGAL::polygon_area_2(poly.begin(), poly.end(), K()));
      CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR;

      if(orientation == expected)
        polys.emplace_back(poly.begin(), poly.end());
      else
        polys.emplace_back(poly.rbegin(), poly.rend());
    }
  };

  if(is_number_of_CC_in_input)
  {
    int ccb_count = 0;
    in >> ccb_count;
    for(int i=0; i<ccb_count && in; ++i)
      read_polygon(i);
  }
  else
  {
    int i = 0;
    while(in)
      read_polygon(i++);
  }

  if(polys.empty())
  {
    std::cerr << "Error: empty input?" << std::endl;
    return false;
  }

  std::cout <<"Polygon with border of size: " << polys[0].size() << std::endl;
  if(polys.size() > 1)
    std::cout << polys.size() - 1 << " hole(s)" << std::endl;

  p = Polygon_with_holes_2(polys[0]);
  for(std::size_t i=0; i<polys.size()-1; ++i)
    p.add_hole(polys[i+1]);

  return true;
}

template <typename K>
bool read_input_polygon(const char* filename,
                        CGAL::Polygon_with_holes_2<K>& p)
{
  std::string ext = CGAL::IO::internal::get_file_extension(filename);
  if(ext == "dat")
  {
    return read_dat_polygon(filename, p);
  }
  else
  {
    std::cerr << "Error: unknown file extension: " << ext << std::endl;
    return false;
  }
}

template <typename K>
bool read_segment_speeds(const char* filename,
                         std::vector<std::vector<typename K::FT> >& weights)
{
  using FT = typename K::FT;

  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Error: Could not read " << filename << std::endl;
    return false;
  }

  std::vector<FT> border_weights;

  std::string line;
  while(getline(in, line))
  {
    if(line.empty())
    {
      weights.push_back(border_weights);
      border_weights.clear();
    }

    std::istringstream iss(line);
    double w;
    if(iss >> w)
      border_weights.push_back(w);
  }

  // in case the last line is not empty
  if(!border_weights.empty())
    weights.push_back(border_weights);

  return true;
}

std::string root_name(const std::string& full_filename)
{
  std::string name = std::string(full_filename);
  name = name.substr(name.find_last_of("/") + 1, name.length() - 1);
  name = name.substr(0, name.find_last_of("."));
  return name;
}

template <typename K>
bool test(const char* poly_filename,
          const char* angles_filename,
          const typename K::FT height,
          const typename K::FT expected_volume)
{
  using FT = typename K::FT;

  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Point_3 = typename K::Point_3;
  using Mesh = CGAL::Surface_mesh<Point_3>;

  std::cout << "\nTest:\n"
            << "\tPolygon: " << poly_filename << "\n"
            << "\tAngles: " << angles_filename << "\n"
            << "\tHeight: " << height << std::endl;

  Polygon_with_holes_2 pwh;
  if(!read_input_polygon(poly_filename, pwh) || pwh.outer_boundary().is_empty())
  {
    std::cerr << "Error: failure during polygon read" << std::endl;
    return false;
  }

  std::vector<std::vector<FT> > angles;
  if(!read_segment_speeds<K>(angles_filename, angles))
  {
    std::cerr << "Error: failure during speed read" << std::endl;
    return false;
  }

  Mesh sm;
  if(pwh.number_of_holes() == 0) // just to test the hole-less version
    CGAL::extrude_skeleton(pwh.outer_boundary(), sm, CGAL::parameters::angles(angles).maximum_height(height));
  else
    CGAL::extrude_skeleton(pwh, sm, CGAL::parameters::angles(angles).maximum_height(height));

//  assert(CGAL::IO::write_polygon_mesh(root_name(poly_filename) + "_extruded_up.off", sm, CGAL::parameters::stream_precision(17)));

  FT volume = PMP::volume(sm);
  const FT rel_eps = 1e-5;

  assert(are_equal(volume, expected_volume, rel_eps, true /*verbose*/));

  // also test with the opposite weight
  clear(sm);
  CGAL::extrude_skeleton(pwh, sm, CGAL::parameters::angles(angles).maximum_height(- height));

//  assert(CGAL::IO::write_polygon_mesh(root_name(poly_filename) + "_extruded_down.off", sm, CGAL::parameters::stream_precision(17)));

  volume = PMP::volume(sm);
  assert(are_equal(volume, expected_volume, rel_eps, true /*verbose*/));

  return true;
}

template <typename K>
void test()
{
  test<K>("data/polygon_000.dat", "data/angles_000.dat",   6, 162.37987499999997);
  test<K>("data/polygon_001.dat", "data/angles_001.dat",   6, 761.76899999999989);
  test<K>("data/polygon_002.dat", "data/angles_002.dat",  22, 15667.658890389464);
  test<K>("data/polygon_003.dat", "data/angles_003.dat",  12, 105.79864291299999);
  test<K>("data/polygon_004.dat", "data/angles_004.dat",  12, 3119.9357499857151);
  test<K>("data/polygon_005.dat", "data/angles_005.dat",  12, 1342.9474791424002);
  test<K>("data/polygon_006.dat", "data/angles_006.dat",  12, 249.41520000000008);
  test<K>("data/polygon_007.dat", "data/angles_007.dat",  12, 7344.8312073148918);
  test<K>("data/polygon_008.dat", "data/angles_008.dat",  12, 7240.6890039677555);
  test<K>("data/polygon_009.dat", "data/angles_009.dat",  10, 3704.0787987580375);
  test<K>("data/polygon_010.dat", "data/angles_010.dat",  40, 29306.453333333335);
  test<K>("data/polygon_011.dat", "data/angles_011.dat",  40, 375866.54633629462);
  test<K>("data/polygon_012.dat", "data/angles_012.dat",  12, 4560.0268861722925);
  test<K>("data/polygon_013.dat", "data/angles_013.dat",  12, 2221.3622594712501);
  test<K>("data/polygon_014.dat", "data/angles_014.dat",  12, 4534.0568515270861);
  test<K>("data/polygon_015.dat", "data/angles_015.dat",  12, 1565.5667825255343);
  test<K>("data/polygon_016.dat", "data/angles_016.dat", 311, 2518611984.6277928);
  test<K>("data/polygon_017.dat", "data/angles_017.dat",  50, 7729166.666666667);
  test<K>("data/polygon_018.dat", "data/angles_018.dat",  50, 354166.66666666663);
  test<K>("data/polygon_019.dat", "data/angles_019.dat", 311, 92570921.033775225);
  test<K>("data/polygon_020.dat", "data/angles_020.dat",  70, 1550161.8131298050);
  test<K>("data/polygon_021.dat", "data/angles_021.dat",  70, 37631800.885042846);
  test<K>("data/polygon_022.dat", "data/angles_022.dat",  20, 7702444.2118858183);
}

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  test<EPICK>();
  test<EPECK>();
  test<EPECK_w_SQRT>();

  std::cout << "OK" << std::endl;

  return EXIT_SUCCESS;
}

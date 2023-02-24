#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/extrude_skeleton.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

namespace PMP = ::CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = K::FT;
using Point_2 = K::Point_2;
using Point_3 = K::Point_3;

using Polygon_2 = CGAL::Polygon_2<K>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

using Straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
using Straight_skeleton_2_ptr = boost::shared_ptr<Straight_skeleton_2>;

using Mesh = CGAL::Surface_mesh<Point_3>;

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

bool read_dat_polygon(const char* filename,
                      Polygon_with_holes_2& p)
{
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

bool read_input_polygon(const char* filename,
                        Polygon_with_holes_2& p)
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

bool read_segment_speeds(const char* filename,
                        std::vector<std::vector<FT> >& weights)
{
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

bool test(const char* poly_filename,
          const char* angles_filename,
          const FT height,
          const FT expected_volume)
{
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
  if(!read_segment_speeds(angles_filename, angles))
  {
    std::cerr << "Error: failure during speed read" << std::endl;
    return false;
  }

  Mesh sm;
  extrude_skeleton(pwh, height, sm, CGAL::parameters::angles(angles));

  const FT volume = PMP::volume(sm);
  const FT rel_eps = 1e-5;

  assert(are_equal(volume, expected_volume, rel_eps, true));

  return true;
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  if(argc != 1)
  {
    Mesh sm;
    CGAL::IO::read_polygon_mesh(argv[1], sm);
    std::cout << PMP::volume(sm) << std::endl;
    return EXIT_SUCCESS;
  }

  test("data/extrusion/polygon_000.dat", "data/extrusion/angles_000.dat",   6, 162.37987499999997);
  test("data/extrusion/polygon_001.dat", "data/extrusion/angles_001.dat",   6, 761.76899999999989);
  test("data/extrusion/polygon_002.dat", "data/extrusion/angles_002.dat",  22, 15667.658890389464);
  test("data/extrusion/polygon_003.dat", "data/extrusion/angles_003.dat",  12, 105.79864291299999);
  test("data/extrusion/polygon_004.dat", "data/extrusion/angles_004.dat",  12, 3119.9357499857151);
  test("data/extrusion/polygon_005.dat", "data/extrusion/angles_005.dat",  12, 1342.9474791424002);
  test("data/extrusion/polygon_006.dat", "data/extrusion/angles_006.dat",  12, 249.41520000000008);
  test("data/extrusion/polygon_007.dat", "data/extrusion/angles_007.dat",  12, 7344.8312073148918);
  test("data/extrusion/polygon_008.dat", "data/extrusion/angles_008.dat",  12, 7240.6890039677555);
  test("data/extrusion/polygon_009.dat", "data/extrusion/angles_009.dat",  10, 3704.0787987580375);
  test("data/extrusion/polygon_010.dat", "data/extrusion/angles_010.dat",  40, 29306.453333333335);
  test("data/extrusion/polygon_011.dat", "data/extrusion/angles_011.dat",  40, 375866.54633629462);
  test("data/extrusion/polygon_012.dat", "data/extrusion/angles_012.dat",  12, 4560.0268861722925);
  test("data/extrusion/polygon_013.dat", "data/extrusion/angles_013.dat",  12, 2221.3622594712501);
  test("data/extrusion/polygon_014.dat", "data/extrusion/angles_014.dat",  12, 4534.0568515270861);
  test("data/extrusion/polygon_015.dat", "data/extrusion/angles_015.dat",  12, 1565.5667825255343);
  test("data/extrusion/polygon_016.dat", "data/extrusion/angles_016.dat", 311, 2518611984.6277928);
  test("data/extrusion/polygon_017.dat", "data/extrusion/angles_017.dat",  50, 7729166.666666667);
  test("data/extrusion/polygon_018.dat", "data/extrusion/angles_018.dat",  50, 354166.66666666663);
  test("data/extrusion/polygon_019.dat", "data/extrusion/angles_019.dat", 311, 92570921.033775225);
  test("data/extrusion/polygon_020.dat", "data/extrusion/angles_020.dat",  70, 1550161.8131298050);
  test("data/extrusion/polygon_021.dat", "data/extrusion/angles_021.dat",  70, 37631800.885042846);
  test("data/extrusion/polygon_022.dat", "data/extrusion/angles_022.dat",  20, 7702444.2118858183);

  std::cout << "OK" << std::endl;

  return EXIT_SUCCESS;
}

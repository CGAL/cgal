#define CGAL_SLS_DEBUG_DRAW

#include <iostream>
#include <iomanip>
#include <string>

// #define CGAL_SLS_PRINT_QUEUE_BEFORE_EACH_POP
// #define CGAL_STRAIGHT_SKELETON_ENABLE_TRACE 100
// #define CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE 10000000
// #define CGAL_STRAIGHT_SKELETON_VALIDITY_ENABLE_TRACE
// #define CGAL_POLYGON_OFFSET_ENABLE_TRACE 10000000

void Straight_skeleton_external_trace(std::string m)
{
  std::cout << std::setprecision(17) << m << std::endl << std::endl ;
}

void Straight_skeleton_traits_external_trace(std::string m)
{
  std::cout << std::setprecision(17) << m << std::endl << std::endl ;
}

#include <CGAL/extrude_skeleton.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/random_polygon_2.h>

#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/draw_surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

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
void test_fixed_with_vertical_combinations()
{
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Point_3 = typename K::Point_3;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Mesh = CGAL::Surface_mesh<Point_3>;

  auto test_combination = [&](const auto& pwh,
                              const std::vector<std::vector<FT> >& angles,
                              const FT max_height)
  {
    Mesh sm;
    bool success = extrude_skeleton(pwh, sm, CGAL::parameters::angles(angles).maximum_height(max_height).verbose(true));
    assert(success);
    if(!success)
    {
      std::cerr << "Error: failed to extrude skeleton" << std::endl;
      std::exit(1);
    }

    CGAL::draw(sm);
    CGAL::IO::write_OFF("last_extrusion.off", sm, CGAL::parameters::stream_precision(17));

    assert(is_closed(sm));
  };

  auto test_all_vertical_combinations = [&](const Polygon_2& poly)
  {
    const std::size_t n = poly.size();
    const FT base_angle = 45;
    const FT vertical_angle = 90;
    const FT height = 2.0; // @todo multiple heights

    // Generate all possible combinations of vertical/non-vertical edges
    for(std::size_t mask = 0; mask < (1u << n); ++mask)
    {
      std::vector<std::vector<FT>> angles;
      std::vector<FT> contour_angles(n);

      // Set angles based on the current mask
      for(std::size_t i = 0; i < n; ++i)
        contour_angles[i] = (mask & (1u << i)) ? vertical_angle : base_angle;

      angles.push_back(contour_angles);

      std::cout << "Testing combination:";
      for(FT angle : contour_angles)
        std::cout << " " << angle;
      std::cout << std::endl;

      test_combination(poly, angles, height);
    }
  };

  std::cout << "\nTesting 'square' combinations:" << std::endl;
  std::vector<Point_2> square_vertices = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
  Polygon_2 square { square_vertices.begin(), square_vertices.end() };
  // test_all_vertical_combinations(square); // @tmp

  std::cout << "\nTesting 'touching squares' combinations:" << std::endl;
  std::vector<Point_2> touching_squares_vertices = {
    { 0, 1 }, { 2, 1 }, { 2, 0 }, { 4, 0 }, { 4, 2 }, { 2, 2 }, { 2, 3 }, { 0, 3 }
  };
  Polygon_2 touching_squares { touching_squares_vertices.begin(), touching_squares_vertices.end() };
  test_all_vertical_combinations(touching_squares);
}

template <typename K>
void test_fixed()
{
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Point_3 = typename K::Point_3;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Mesh = CGAL::Surface_mesh<Point_3>;

  auto test_fixed_polygon = [&](const auto& pwh,
                                const std::vector<std::vector<FT> >& angles,
                                const FT max_height,
                                const FT expected_vol)
  {
    Mesh sm;
    bool success = extrude_skeleton(pwh, sm, CGAL::parameters::angles(angles).verbose(true).maximum_height(max_height));
    assert(success);
    if(!success)
    {
      std::cerr << "Error: failed to extrude skeleton" << std::endl;
      std::exit(1);
    }

    CGAL::draw(sm);
    CGAL::IO::write_OFF("last_extrusion.off", sm, CGAL::parameters::stream_precision(17));

    FT actual_vol = CGAL::Polygon_mesh_processing::volume(sm);
    std::cout << "Expected volume: " << expected_vol << ", actual is: " << actual_vol << std::endl;
    if (CGAL::abs(expected_vol - actual_vol) > 1e-6)
    {
      assert(false);
      std::exit(1);
    }
  };

  // unit square
  std::vector<Point_2> square_vertices = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
  Polygon_2 square { square_vertices.begin(), square_vertices.end() };

  // unit square, 45° inward extrusion
  const std::vector<std::vector<FT> > square_angles_45 = { { 45, 45, 45, 45 } };
  test_fixed_polygon(square, square_angles_45,  1, FT(1)/6);
  test_fixed_polygon(square, square_angles_45, -1, FT(1)/6);

  // unit square, 80° inward extrusion
  const std::vector<std::vector<FT> > square_angles_80 = { { 80, 80, 80, 80 } };
  test_fixed_polygon(square, square_angles_80,  1000000, 0.9452136366);
  test_fixed_polygon(square, square_angles_80, -1000000, 0.9452136366);
  test_fixed_polygon(square, square_angles_80,        1, 0.6888009778);
  test_fixed_polygon(square, square_angles_80,       -1, 0.6888009778);

  // unit square, 10° outward extrusion
  // the volume of a truncated pyramid is: 1/3 × h × (a^2 + b^2 + ab)
  // and the extruded square has side length: b = a + 2 * (h * tan(alpha - 90°))
  const std::vector<std::vector<FT> > square_angles_m10 = { { 170, 170, 170, 170 } };
  test_fixed_polygon(square, square_angles_m10,   1, 55.2271469426);
  test_fixed_polygon(square, square_angles_m10,  -1, 55.2271469426);
  test_fixed_polygon(square, square_angles_m10,  10, 44028.8396672);
  test_fixed_polygon(square, square_angles_m10, -10, 44028.8396672);

  const std::vector<std::vector<FT> > square_angles_m45 = { { 135, 135, 135, 135 } };
  test_fixed_polygon(square, square_angles_m45,   1, FT(13)/3);
  test_fixed_polygon(square, square_angles_m45,  -1, FT(13)/3);
  test_fixed_polygon(square, square_angles_m45,  10, FT(4630)/3);
  test_fixed_polygon(square, square_angles_m45, -10, FT(4630)/3);

  // some vertical
  const std::vector<std::vector<FT> > square_semi_vertical_1 = { { 90, 45, 45, 45 } };
  test_fixed_polygon(square, square_semi_vertical_1,   1,  FT(1)/6 + FT(1)/24);
  test_fixed_polygon(square, square_semi_vertical_1,  -1,  FT(1)/6 + FT(1)/24);

  const std::vector<std::vector<FT> > square_semi_vertical_2 = { { 90, 90, 45, 45 } };
  test_fixed_polygon(square, square_semi_vertical_2,   1,  FT(1)/3);
  test_fixed_polygon(square, square_semi_vertical_2,  -1,  FT(1)/3);

  const std::vector<std::vector<FT> > square_semi_vertical_2b = { { 90, 45, 90, 45 } };
  test_fixed_polygon(square, square_semi_vertical_2b,   1, 0.25);
  test_fixed_polygon(square, square_semi_vertical_2b,  -1, 0.25);

  const std::vector<std::vector<FT> > square_semi_vertical_3 = { { 90, 90, 90, 45 } };
  test_fixed_polygon(square, square_semi_vertical_3,   1, FT(1)/2);
  test_fixed_polygon(square, square_semi_vertical_3,  -1, FT(1)/2);

  // unit square, 10° inward extrusion
  const std::vector<std::vector<FT> > square_angles_10 = { { 10, 10, 10, 10 } };
  test_fixed_polygon(square, square_angles_10,  1, 0.0293878301); // 0.5 * tan(alpha) / 3
  test_fixed_polygon(square, square_angles_10, -1, 0.0293878301);

  // all vertical
  const std::vector<std::vector<FT> > square_angles_90 = { { 90, 90, 90, 90 } };
  test_fixed_polygon(square, square_angles_90,   1,  1);
  test_fixed_polygon(square, square_angles_90,  -1,  1);
  test_fixed_polygon(square, square_angles_90,  10, 10);
  test_fixed_polygon(square, square_angles_90, -10, 10);

  // larger square, with holes
  square_vertices = { { 0, 0 }, { 10, 0 }, { 10, 10 }, { 0, 10 } };
  Polygon_2 large_square { square_vertices.begin(), square_vertices.end() };

  Polygon_with_holes_2 pwh(large_square);
  std::vector<Point_2> hole_centers = { { 2.5, 2.5 }, { 7.5, 2.5 }, { 7.5, 7.5 }, { 2.5, 7.5 } };
  for(const Point_2& hole_center : hole_centers)
  {
    std::vector<Point_2> hole_pts = { { hole_center.x() - 0.5, hole_center.y() - 0.5 },
                                      { hole_center.x() - 0.5, hole_center.y() + 0.5 },
                                      { hole_center.x() + 0.5, hole_center.y() + 0.5 },
                                      { hole_center.x() + 0.5, hole_center.y() - 0.5 } };
    Polygon_2 hole(hole_pts.begin(), hole_pts.end());
    pwh.add_hole(hole);
  }

  // all vertical
  const std::vector<std::vector<FT> > large_square_angles_90 = { { 90, 90, 90, 90 }, { 90, 90, 90, 90 }, { 90, 90, 90, 90 }, { 90, 90, 90, 90 }, { 90, 90, 90, 90 } };
  test_fixed_polygon(pwh, large_square_angles_90,   1, 100 - 4 * 1);
  test_fixed_polygon(pwh, large_square_angles_90,  -1, 100 - 4 * 1);
  test_fixed_polygon(pwh, large_square_angles_90,  10, 1000 - 4 * 10);
  test_fixed_polygon(pwh, large_square_angles_90, -10, 1000 - 4 * 10);

  // outer vertical, holes outward
  const std::vector<std::vector<FT> > large_square_angles_outer_90 = { { 90, 90, 90, 90 }, { 135, 135, 135, 135 }, { 135, 135, 135, 135 }, { 135, 135, 135, 135 }, { 135, 135, 135, 135 } };
  test_fixed_polygon(pwh, large_square_angles_outer_90,   1, 100 - 4 * FT(1)/6);
  test_fixed_polygon(pwh, large_square_angles_outer_90,  -1, 100 - 4 * FT(1)/6);
  test_fixed_polygon(pwh, large_square_angles_outer_90,  10, 1000 - 4 * FT(1)/6);
  test_fixed_polygon(pwh, large_square_angles_outer_90, -10, 1000 - 4 * FT(1)/6);

  // outer vertical, holes inwards
  const std::vector<std::vector<FT> > large_square_angles_outer_90b = { { 90, 90, 90, 90 }, { 45, 45, 45, 45 }, { 45, 45, 45, 45 }, { 45, 45, 45, 45 }, { 45, 45, 45, 45 } };
  test_fixed_polygon(pwh, large_square_angles_outer_90b,   1, 100 - 4 * FT(13) / 3); // see 45° above
  test_fixed_polygon(pwh, large_square_angles_outer_90b,  -1, 100 - 4 * FT(13) / 3);
  test_fixed_polygon(pwh, large_square_angles_outer_90b,  10, 200 - (4 * FT(62) / 3)); // h = 2 since 45°
  test_fixed_polygon(pwh, large_square_angles_outer_90b, -10, 200 - (4 * FT(62) / 3));

  // outer inwards, holes vertical
  const std::vector<std::vector<FT> > large_square_angles_some_90 = { { 10, 10, 10, 10 }, { 90, 90, 90, 90 }, { 90, 90, 90, 90 }, { 90, 90, 90, 90 }, { 90, 90, 90, 90 } };
  test_fixed_polygon(pwh, large_square_angles_some_90,   0.1, 8.5086282194);
  test_fixed_polygon(pwh, large_square_angles_some_90,  -0.1, 8.5086282194);
  test_fixed_polygon(pwh, large_square_angles_some_90,    10, 27.742111630);
  test_fixed_polygon(pwh, large_square_angles_some_90,   -10, 27.742111630);

  // outer outwardcs, holes vertical
  const std::vector<std::vector<FT> > large_square_angles_some_90b = { { 170, 170, 170, 170 }, { 90, 90, 90, 90 }, { 90, 90, 90, 90 }, { 90, 90, 90, 90 }, { 90, 90, 90, 90 } };
  test_fixed_polygon(pwh, large_square_angles_some_90b,   1, 256.310219695 - 4 * 1);
  test_fixed_polygon(pwh, large_square_angles_some_90b,  -1, 256.310219695 - 4 * 1);
  test_fixed_polygon(pwh, large_square_angles_some_90b,  10, 55227.1469423 - 4 * 10);
  test_fixed_polygon(pwh, large_square_angles_some_90b, -10, 55227.1469423 - 4 * 10);

  // everyone outwards
  const std::vector<std::vector<FT> > large_square_angles_all_170 = { { 170, 170, 170, 170 }, { 170, 170, 170, 170 }, { 170, 170, 170, 170 }, { 170, 170, 170, 170 }, { 170, 170, 170, 170 } };
  test_fixed_polygon(pwh, large_square_angles_all_170,   1, 256.310219695 - 4 * 0.0293878301); // see 10° above
  test_fixed_polygon(pwh, large_square_angles_all_170,  -1, 256.310219695 - 4 * 0.0293878301);
  test_fixed_polygon(pwh, large_square_angles_all_170,  10, 55227.1469423 - 4 * 0.0293878301);
  test_fixed_polygon(pwh, large_square_angles_all_170, -10, 55227.1469423 - 4 * 0.0293878301);

  // everyone inwards
  const std::vector<std::vector<FT> > large_square_angles_all_45 = { { 45, 45, 45, 45 }, { 45, 45, 45, 45 }, { 45, 45, 45, 45 }, { 45, 45, 45, 45 }, { 45, 45, 45, 45 } };
  test_fixed_polygon(pwh, large_square_angles_all_45,   1, 64); // big truncated pyramid - 4 small truncated pyramids
  test_fixed_polygon(pwh, large_square_angles_all_45,  -1, 64);
  test_fixed_polygon(pwh, large_square_angles_all_45,  10, FT(232)/3);
  test_fixed_polygon(pwh, large_square_angles_all_45, -10, FT(232)/3);
}

template <typename K>
bool test_dat(const char* poly_filename,
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

//  CGAL::IO::write_polygon_mesh(root_name(poly_filename) + "_extruded_up.off", sm, CGAL::parameters::stream_precision(17));

  CGAL::draw(sm);

  FT volume = PMP::volume(sm);
  const FT rel_eps = 1e-5;

  assert(are_equal(volume, expected_volume, rel_eps, true /*verbose*/));

  // also test with the opposite weight
  clear(sm);
  CGAL::extrude_skeleton(pwh, sm, CGAL::parameters::angles(angles).maximum_height(- height));

//  CGAL::IO::write_polygon_mesh(root_name(poly_filename) + "_extruded_down.off", sm, CGAL::parameters::stream_precision(17));

  volume = PMP::volume(sm);
  assert(are_equal(volume, expected_volume, rel_eps, true /*verbose*/));

  return true;
}


template <typename K>
void test_random(CGAL::Random& rnd)
{
  std::cout << "\n ==== Test with Kernel: " << typeid(K).name() << " ====" << std::endl;

  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Point_3 = typename K::Point_3;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Mesh = CGAL::Surface_mesh<Point_3>;

  // ------------------------ test some random polygons ------------------------
  auto generate_random_polygon = [&](CGAL::Random& rnd) -> Polygon_2
  {
    typedef CGAL::Random_points_in_square_2<Point_2> Point_generator;

    Polygon_2 poly;
    CGAL::random_polygon_2(10, std::back_inserter(poly), Point_generator(0.25, rnd));
    return poly;
  };

  auto test_polygon = [&](const auto& pwh)
  {
    std::cout << "== Test Polygon ==" << std::endl;
    for(const auto& p : pwh.outer_boundary())
      std::cout << "  " << p << std::endl;

    // test all slopes vertical
    {
      for(int i=0; i<3; ++i)
      {
        const FT max_height = rnd.get_double(-1, 1);
        if(max_height == 0)
          continue;

        std::cout << "max_height = " << max_height << std::endl;

        std::vector<std::vector<FT> > weights;
        weights.push_back(std::vector<FT>(pwh.outer_boundary().size(), 0));
        for(auto hit=pwh.holes_begin(); hit!=pwh.holes_end(); ++hit)
          weights.push_back(std::vector<FT>(hit->size(), 0));

        Mesh sm;
        bool success = extrude_skeleton(pwh, sm, CGAL::parameters::weights(weights).verbose(true).maximum_height(max_height));
        assert(success);
        if(!success)
        {
          std::cerr << "Error: failed to extrude skeleton" << std::endl;
          assert(false);
          std::exit(1);
        }

        CGAL::draw(sm);
        CGAL::IO::write_OFF("last_extrusion.off", sm, CGAL::parameters::stream_precision(17));

        // check the volume
        FT expected_vol = CGAL::abs(pwh.outer_boundary().area() * max_height);
        FT actual_vol = CGAL::Polygon_mesh_processing::volume(sm);
        if (CGAL::abs(expected_vol - actual_vol) > 1e-6)
        {
          std::cerr << "Error: expected volume " << expected_vol << " but got " << actual_vol << std::endl;
          assert(false);
          std::exit(1);
        }

        // for each point at z=0, there should be the same point at z=max_height
        std::map<std::pair<FT, FT>, std::list<Point_3> > zs;
        for (auto v : vertices(sm))
        {
          Point_3 p = sm.point(v);
          if (p.z() != 0 && p.z() != max_height)
          {
            std::cerr << "Error: point at z=" << p.z() << " instead of 0 or " << max_height << std::endl;
            assert(false);
            std::exit(1);
          }

          zs[std::make_pair(p.x(), p.y())].push_back(p);
        }

        for (const auto& e : zs)
        {
          if (e.second.size() != 2)
          {
            std::cerr << "Error: point at (" << e.first.first << ", " << e.first.second << ") has " << e.second.size() << " points" << std::endl;
            assert(false);
            std::exit(1);
          }

          if (e.second.front().z() == e.second.back().z())
          {
            std::cerr << "Error: point at (" << e.first.first << ", " << e.first.second << ") has two points at the same height" << std::endl;
            assert(false);
            std::exit(1);
          }
        }
      }
    }

    // test some slopes vertical
    {
      for(int i=0; i<3; ++i)
      {
        const FT max_height = rnd.get_double(-1, 1);
        if(max_height == 0)
          continue;

        std::vector<std::vector<FT> > weights;
        weights.push_back(std::vector<FT>(pwh.outer_boundary().size(), 0));
        for(auto hit=pwh.holes_begin(); hit!=pwh.holes_end(); ++hit)
          weights.push_back(std::vector<FT>(hit->size(), 0));

        for(auto& contour_weights : weights)
        {
          for(FT& weight : contour_weights)
          {
            weight = rnd.get_double(0.05, 0.5);
            if (!rnd.get_int(0, 3)) { // a third to zero
              weight = 0;
            }
          }
        }

        // randomly switch inwards/outwards
        bool outwards = false;
        if(rnd.get_int(0, 2))
        {
          outwards = true;
          for(std::vector<FT>& ws : weights)
            for(FT& w : ws)
              w = -w;
        }

        Mesh sm;
        bool success = extrude_skeleton(pwh, sm, CGAL::parameters::weights(weights)
                                                                  .verbose(true)
                                                                  .maximum_height(max_height));
        if(!success)
        {
          std::cerr << "Error: failed to extrude skeleton" << std::endl;
          assert(false);
          std::exit(1);
        }

        CGAL::draw(sm);
        CGAL::IO::write_OFF("last_extrusion.off", sm, CGAL::parameters::stream_precision(17));

        if(!outwards)
        {
          success = extrude_skeleton(pwh, sm, CGAL::parameters::weights(weights).verbose(true));
          if(!success)
          {
            std::cerr << "Error: failed to extrude skeleton" << std::endl;
            assert(false);
            std::exit(1);
          }

          CGAL::draw(sm);
          CGAL::IO::write_OFF("last_extrusion.off", sm, CGAL::parameters::stream_precision(17));
        }
      }
    }
  };

  // Random simple polygon (no holes)
  for (int i=0; i<10; ++i) {
    Polygon_2 poly = generate_random_polygon(rnd);
    Polygon_with_holes_2 pwh(poly);
    test_polygon(pwh);
  }

  // Random simple polygons with holes

}

template <typename K>
void test_dats()
{
  test_dat<K>("data/polygon_000.dat", "data/angles_000.dat",   6, 162.37987499999997);
  test_dat<K>("data/polygon_001.dat", "data/angles_001.dat",   6, 761.76899999999989);
  test_dat<K>("data/polygon_002.dat", "data/angles_002.dat",  22, 15667.658890389464);
  test_dat<K>("data/polygon_003.dat", "data/angles_003.dat",  12, 105.79864291299999);
  test_dat<K>("data/polygon_004.dat", "data/angles_004.dat",  12, 3119.9357499857151);
  test_dat<K>("data/polygon_005.dat", "data/angles_005.dat",  12, 1342.9474791424002);
  test_dat<K>("data/polygon_006.dat", "data/angles_006.dat",  12, 249.41520000000008);
  test_dat<K>("data/polygon_007.dat", "data/angles_007.dat",  12, 7344.8312073148918);
  test_dat<K>("data/polygon_008.dat", "data/angles_008.dat",  12, 7240.6890039677555);
  test_dat<K>("data/polygon_009.dat", "data/angles_009.dat",  10, 3704.0787987580375);
  test_dat<K>("data/polygon_010.dat", "data/angles_010.dat",  40, 29306.453333333335);
  test_dat<K>("data/polygon_011.dat", "data/angles_011.dat",  40, 375866.54633629462);
  test_dat<K>("data/polygon_012.dat", "data/angles_012.dat",  12, 4560.0268861722925);
  test_dat<K>("data/polygon_013.dat", "data/angles_013.dat",  12, 2221.3622594712501);
  test_dat<K>("data/polygon_014.dat", "data/angles_014.dat",  12, 4534.0568515270861);
  test_dat<K>("data/polygon_015.dat", "data/angles_015.dat",  12, 1565.5667825255343);
  test_dat<K>("data/polygon_016.dat", "data/angles_016.dat", 311, 2518611984.6277928);
  test_dat<K>("data/polygon_017.dat", "data/angles_017.dat",  50, 7729166.666666667);
  test_dat<K>("data/polygon_018.dat", "data/angles_018.dat",  50, 354166.66666666663);
  test_dat<K>("data/polygon_019.dat", "data/angles_019.dat", 311, 92570921.033775225);
  test_dat<K>("data/polygon_020.dat", "data/angles_020.dat",  70, 1550161.8131298050);
  test_dat<K>("data/polygon_021.dat", "data/angles_021.dat",  70, 37631800.885042846);
  test_dat<K>("data/polygon_022.dat", "data/angles_022.dat",  20, 7702444.2118858183);
}


template <typename K>
void test(CGAL::Random& rnd)
{
  test_fixed_with_vertical_combinations<K>();
  // test_fixed<K>();
  // test_dats<K>();
  // test_random<K>(rnd);
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  int seed = (argc > 0) ? std::atoi(argv[0]) : 0; // @fixme std::time(nullptr)

  CGAL::Random rnd(seed);
  std::cout << "Seed is " << rnd.get_seed() << std::endl;

  test<EPICK>(rnd);
  test<EPECK>(rnd);
  test<EPECK_w_SQRT>(rnd);

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}

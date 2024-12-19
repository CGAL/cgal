
#include "test_util.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/internal/implicit_shapes_helper.h>

#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/Random.h>

#include <cassert>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#define CGAL_TESTSUITE_ISOSURFACING_OUTPUT

namespace IS = CGAL::Isosurfacing;

template <typename Grid>
bool check_closed_not_empty(const Grid& grid, const IS::Interpolated_discrete_values_3<Grid>& values, const typename Grid::Geom_traits::FT iso = 0)
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  FT boundary_min = std::numeric_limits<FT>::max();
  FT total_min = std::numeric_limits<FT>::max();
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        total_min = (std::min)(total_min, values(x, y, z));
        if (x == 0 || y == 0 || z == 0 || x == grid.xdim() - 1 || y == grid.ydim() - 1 || z == grid.zdim() - 1)
          boundary_min = (std::min)(boundary_min, values(x, y, z));
      }
    }
  }

  // empty cell
  if (total_min > iso)
  {
    return false;
  }

  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        const bool is_boundary = x == 0 || y == 0 || z == 0 || x == grid.xdim() - 1 || y == grid.ydim() - 1 || z == grid.zdim() - 1;

        // closed field
        if (!is_boundary && boundary_min - values(x, y, z) < 1e-6)
          return false;
      }
    }
  }
  return true;
}

template <typename Grid>
bool check_iso_vertices(const Grid& grid, const IS::Interpolated_discrete_values_3<Grid>& values, const typename Grid::Geom_traits::FT iso = 0) {
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (std::abs(values(x, y, z)) < iso + 1e-6)
          return true;
      }
    }
  }
  return false;
}

template <typename Grid>
bool check_iso_edges(const Grid& grid, const IS::Interpolated_discrete_values_3<Grid>& values, const typename Grid::Geom_traits::FT iso = 0) {
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (x < grid.xdim() - 1 && std::abs(values(x, y, z)) < iso + 1e-6 && std::abs(values(x + 1, y, z)) < iso + 1e-6)
          return true;
        if (y < grid.ydim() - 1 && std::abs(values(x, y, z)) < iso + 1e-6 && std::abs(values(x, y + 1, z)) < iso + 1e-6)
          return true;
        if (z < grid.zdim() - 1 && std::abs(values(x, y, z)) < iso + 1e-6 && std::abs(values(x, y, z + 1)) < iso + 1e-6)
          return true;
      }
    }
  }
  return false;
}

template <typename Grid>
void to_grid(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values, const std::array<typename Grid::Geom_traits::FT, 8>& center_values) {
  assert(grid.xdim() == 2);
  assert(grid.ydim() == 2);
  assert(grid.zdim() == 2);

  values(0, 0, 0) = center_values[0];
  values(0, 0, 1) = center_values[1];
  values(0, 1, 0) = center_values[2];
  values(0, 1, 1) = center_values[3];
  values(1, 0, 0) = center_values[4];
  values(1, 0, 1) = center_values[5];
  values(1, 1, 0) = center_values[6];
  values(1, 1, 1) = center_values[7];
}

template <typename Grid>
void embed_in_positive_grid(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values, const std::array<typename Grid::Geom_traits::FT, 8>& center_values) {
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  assert(grid.xdim() == 4);
  assert(grid.ydim() == 4);
  assert(grid.zdim() == 4);

  CGAL::Random rand;

  FT min = (std::numeric_limits<FT>::max)();
  FT max = std::numeric_limits<FT>::lowest();
  for (const FT v : center_values)
  {
    min = (std::min)(min, v);
    max = (std::max)(max, v);
  }
  max += 1e-3;

  for(std::size_t x=0; x<grid.xdim(); ++x)
    for(std::size_t y=0; y<grid.ydim(); ++y)
      for(std::size_t z=0; z<grid.zdim(); ++z)
        values(x, y, z) = rand.uniform_real<FT>(max, 2 * max - min);

  values(1, 1, 1) = center_values[0];
  values(1, 1, 2) = center_values[1];
  values(1, 2, 1) = center_values[2];
  values(1, 2, 2) = center_values[3];
  values(2, 1, 1) = center_values[4];
  values(2, 1, 2) = center_values[5];
  values(2, 2, 1) = center_values[6];
  values(2, 2, 2) = center_values[7];
}

enum class AmbiguousCase {
  CONTOUR_6_VTS,
  CONTOUR_7_VTS_CONFIG_A,
  CONTOUR_7_VTS_CONFIG_B,
  CONTOUR_8_VTS_CONFIG_A,
  CONTOUR_8_VTS_CONFIG_B,
  CONTOUR_9_VTS_CONFIG_A,
  MC_4_TUNNEL,
  MC_7_TUNNEL,
  MC_10_TUNNEL,
  MC_12_TUNNEL,
  MC_13_TUNNEL,
  MC_13_12_VTS_CONFIG_A,
  MC_13_12_VTS_CONFIG_B,
  MC_13_SIMPLE
};

template <typename FT>
FT generate_predefined_inner_ambiguity(const AmbiguousCase topology_case, std::array<FT, 8>& values) {
  FT iso;

  switch (topology_case) {
    case AmbiguousCase::CONTOUR_6_VTS:
      values = {-0.960492426392903, 0.793207329559554, 0.916735525189067, -0.422761282626275, -0.934993247757551, -0.850129305868777, -0.0367116785741896, -0.656740699156587};
      iso = 0.0387;
      break;
    case AmbiguousCase::CONTOUR_7_VTS_CONFIG_A:
      values = {10.2967816247556, 9.45145192686147, 9.54753711271687, 10.6482067822841, 9.81494966341055, 9.31168538578250, 9.80950580411527, 10.7451536262220};
      iso = 9.8588;
      break;
    case AmbiguousCase::CONTOUR_7_VTS_CONFIG_B:
      values = {9.9998593195995547, 9.9993381282115549, 9.9979160205452544, 9.9986053863704142, 9.9999374908631235, 9.999424800002032, 9.9983922749132219, 9.999579324965488};
      iso = 9.9994608191478135;
      break;
    case AmbiguousCase::CONTOUR_8_VTS_CONFIG_A:
      values = {0.454797708726920, 0.801330575352402, 0.649991492712356, -0.973974554763863,  -0.134171007607162,  -0.0844698148589140, -0.826313795402046, 0.433391503783462};
      iso = 0.0097;
      break;
    case AmbiguousCase::CONTOUR_8_VTS_CONFIG_B:
      values = {9.9985934885536665, 9.9998695572230147, 9.9999045831713928, 9.999316745478131, 9.9986117521866866, 9.9998754368055813, 9.9999031760062458, 9.9992041920402936};
      iso = 9.99946;
      break;
    case AmbiguousCase::CONTOUR_9_VTS_CONFIG_A:
      values = {-15.6504952739285, 2.90290077342601, 24.5454566157887, -24.5274127623786, 21.6741877710053, -4.49696327433901, -19.7891575872492, -15.5588482753161};
      iso = -1.8061;
      break;
    case AmbiguousCase::MC_4_TUNNEL:
      values = {-7.70146936482581, -3.21868369245987, -5.44023748418735, 15.6051950593180, 12.7611835388515, -4.46952393442309, -11.7240576326183, -9.23038948829007};
      iso = -1.7660;
      break;
    case AmbiguousCase::MC_7_TUNNEL:
      values = {-3.42744283804455, 0.621278122151001, 4.48110777981235, -1.95551129669134, 2.30448107596369, -1.04182240925489, -3.51087814405650, -6.44976786808517};
      iso = -0.6;
      break;
    case AmbiguousCase::MC_10_TUNNEL:
      values = {-0.100000000000000, -6.11000000000000, 2, 10.2000000000000, 10.8000000000000, 1.80000000000000, -8.20000000000000, -0.180000000000000};
      iso = 1.10918;
      break;
    case AmbiguousCase::MC_12_TUNNEL:
      values = {-3.37811990337124, 0.473258332744286, 2.54344310345736, 7.87658724379480, 4.38700713005133, -1.49950251870885, -4.21025867362045, -1.00233824192217};
      iso = 0.0708;
      break;
    case AmbiguousCase::MC_13_TUNNEL:
      values = {2.74742516087490, -3.39187542578189, -12.5297639669456, 0.431517989649243, -6.92460546400188, 2.52228314017858, 14.6950568276448, -10.0732624062474};
      iso = -1.3064;
      break;
    case AmbiguousCase::MC_13_12_VTS_CONFIG_A:
      values = {0.546912886195662, -0.421103532406922, -0.643375084081520, 0.855507421818445, -0.260686312588506, 0.206413666735986, 0.237274227130530, -0.183297728364877};
      iso = 0.0293;
      break;
    case AmbiguousCase::MC_13_12_VTS_CONFIG_B:
      values = {1069, 843, 950, 1133, 958, 1029, 1198, 946};
      iso = 1007.4;
      break;
    case AmbiguousCase::MC_13_SIMPLE:
      values = {0.520482995461163, -0.839814387388296, -0.467491517013617, 0.937814095887345, -0.825777099007084, 0.506695544835103, 0.345318915961394, -0.861107217966913};
      iso = 0.0293;
      break;
    default:
      assert(false);
  }
  return iso;
}

enum class SingularCase {
  SIMPLE,
  DOUBLE,
  TRIPLE,
  TUNNEL_OPEN,
  TUNNEL,
  TUNNEL_DOUBLE
};

template <typename FT>
FT generate_predefined_singular(const SingularCase topology_case, std::array<FT, 8>& values) {
  FT iso;
  bool invert;

  switch (topology_case) {
    case SingularCase::SIMPLE:
      values = {-1, 1, -1, -1, -1, -1, -1, 1};
      iso = 0;
      invert = false;
      break;
    case SingularCase::DOUBLE:
      values = {-1, -1, -0.5, 0.5, -2, 1, 0.5, -0.5};
      iso = 0;
      invert = false;
      break;
    case SingularCase::TRIPLE:
      values = {-1, 1, -1, -1, 1, -1, -1, 1};
      iso = 0;
      invert = false;
      break;
    case SingularCase::TUNNEL_OPEN:
      values = {-0.5, 0.5, -0.4, 0.46, 0.5, -0.5, 0.71, -0.92};
      iso = 0;
      invert = false;
      break;
    case SingularCase::TUNNEL:
      values = {-0.5, 0.5, -0.33, 0.46, 0.5, -0.5, 0.71, -0.92};
      iso = 0;
      invert = false;
      break;
    case SingularCase::TUNNEL_DOUBLE:
      values = {-1.2, -0.6, 1.2, 0.9, 9.7, 0.6, -9.7, -0.9};
      iso = 0;
      invert = false;
      break;
    default:
      assert(false);
  }

  for (FT& v : values)
  {
    if (invert)
      v = 2 * iso - v;
  }

  return iso;
}

enum class IsoEdgeCase {
  SIMPLE,
  SEPARATED,
  SINGULAR,
  DOUBLE_SINGULAR_SIMPLE,
  DOUBLE_SINGULAR,
  DOUBLE_SINGULAR_SPLIT,
  WEDGE,
  WEDGE_SINGULAR,
  EDGE_TUNNEL
};

template <typename FT>
FT generate_predefined_iso_edge(const IsoEdgeCase topology_case, std::array<FT, 8>& values) {
  FT iso;
  bool invert;

  switch (topology_case) {
    case IsoEdgeCase::SIMPLE:
      values = {0, 0, -0.4, -0.3, 0.7, 0.5, -0.71, -0.92};
      iso = 0;
      invert = false;
      break;
    case IsoEdgeCase::SEPARATED:
      values = {0, 0, 0.4, 0.3, 0.7, 0.5, -0.71, -0.92};
      iso = 0;
      invert = false;
      break;
    case IsoEdgeCase::SINGULAR:
      values = {0, 0, 0.4, 0.3, -0.7, 0.5, -0.71, -0.92};
      iso = 0;
      invert = false;
      break;
    case IsoEdgeCase::DOUBLE_SINGULAR_SIMPLE:
      values = {0, 0, -0.8, 0.3, -0.7, 0.5, -0.2, 0.1};
      iso = 0;
      invert = false;
      break;
    case IsoEdgeCase::DOUBLE_SINGULAR:
      values = {0, 0, -0.2, 0.3, -0.7, 0.5, -0.2, -0.92};
      iso = 0;
      invert = false;
      break;
    case IsoEdgeCase::DOUBLE_SINGULAR_SPLIT:
      values = {0, 0, 0.2, -0.3, -0.7, 0.5, -0.2, -0.92};
      iso = 0;
      invert = false;
      break;
    case IsoEdgeCase::WEDGE:
      values = {0, 0, -0.4, -0.3, 0, 0.5, -0.2, -0.92};
      iso = 0;
      invert = false;
      break;
    case IsoEdgeCase::WEDGE_SINGULAR:
      values = {0, 0, 0.2, 0.3, 0, 0.5, -0.2, -0.92};
      iso = 0;
      invert = false;
      break;
    case IsoEdgeCase::EDGE_TUNNEL:
      values = {0.3, 0.4, 0.7, -0.4, -0.7, 0, 0.4, 0};
      iso = 0;
      invert = false;
      break;
    default:
      assert(false);
  }

  for (FT& v : values)
  {
    if (invert)
      v = 2 * iso - v;
  }

  return iso;
}

enum class PlaneCase {
  SIMPLE,
  INTERSECTING,
  PLANE_TUBES,
  CROSS,
};

template <typename FT>
FT generate_predefined_plane(const PlaneCase topology_case, std::array<FT, 8>& values) {
  FT iso;
  bool invert;

  switch (topology_case) {
    case PlaneCase::SIMPLE:
      values = {0, 0, 0, 0, 0.7, 0.5, 0.71, 0.92};
      iso = 0;
      invert = false;
      break;
    case PlaneCase::INTERSECTING:
      values = {-12, -91, 12, 91, 97, 9, -97, -9};
      iso = 0;
      invert = false;
      break;
    case PlaneCase::PLANE_TUBES:
      values = {-0.5, 0, 0.4, 0, 0.5, 0, -0.5, 0};
      iso = 0;
      invert = false;
      break;
    case PlaneCase::CROSS:
      values = {1, -1, -1, 1, -1, 1, 1, -1};
      iso = 0;
      invert = false;
      break;
    default:
      assert(false);
  }

  for (FT& v : values)
  {
    if (invert)
      v = 2 * iso - v;
  }

  return iso;
}

template<typename KERNEL>
void compare_tmc_mc_trilinear(const std::array<typename KERNEL::FT, 8>& case_values, typename KERNEL::FT iso)
{
  using K = KERNEL;
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Interpolated_discrete_values_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;
  using Point = typename K::Point_3;
  using Mesh = CGAL::Surface_mesh<Point>;
  using Point_range = std::vector<Point>;
  using Triangle_range = std::vector<std::vector<std::size_t> >;

  Grid grid { Point{-1., -1., -1.}, Point{1., 1., 1.}, std::array<std::size_t, 3>{2, 2, 2} };
  Values values { grid };
  Domain domain { grid, values };

  write_debug_grid(domain, "debug_grid_2x2.off");

  to_grid(grid, values, case_values);

  {
    Grid grid_high_res { Point{-1., -1., -1.}, Point{1., 1., 1.}, std::array<std::size_t, 3>{151, 151, 151} };
    IS::Value_function_3<Grid> values_high_res { values, grid_high_res };
    IS::Marching_cubes_domain_3<Grid, IS::Value_function_3<Grid>> domain_high_res { grid_high_res, values_high_res };

    Point_range points_high_res;
    Triangle_range triangles_high_res;
    IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain_high_res, iso, points_high_res, triangles_high_res, CGAL::parameters::use_topologically_correct_marching_cubes(true));

    #ifdef CGAL_TESTSUITE_ISOSURFACING_OUTPUT
    CGAL::IO::write_polygon_soup("trilinear.off", points_high_res, triangles_high_res);
    #endif

    assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles_high_res));
    Mesh m;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points_high_res, triangles_high_res, m);

    std::cout << "components: " << connected_components(m) << std::endl;
    std::cout << "euler: " << euler_characteristic(m) << std::endl;
    std::cout << "boundaries: " << boundary_components(m) << std::endl;
  }

  {
    Point_range points_mc;
    Triangle_range triangles_mc;
    IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, iso, points_mc, triangles_mc, CGAL::parameters::use_topologically_correct_marching_cubes(false));

    #ifdef CGAL_TESTSUITE_ISOSURFACING_OUTPUT
    CGAL::IO::write_polygon_soup("mc.off", points_mc, triangles_mc);
    #endif
  }

  Point_range points;
  Triangle_range triangles;
  IS::marching_cubes<CGAL::Sequential_tag>(domain, iso, points, triangles, CGAL::parameters::use_topologically_correct_marching_cubes(true));

  #ifdef CGAL_TESTSUITE_ISOSURFACING_OUTPUT
  CGAL::IO::write_polygon_soup("tmc.off", points, triangles);
  #endif
}

template<typename KERNEL>
void assert_tmc(const std::array<typename KERNEL::FT, 8>& case_values, typename KERNEL::FT iso, std::size_t components, std::size_t euler, std::size_t boundaries)
{
  using K = KERNEL;
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Interpolated_discrete_values_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;
  using Point = typename K::Point_3;
  using Mesh = CGAL::Surface_mesh<Point>;
  using Point_range = std::vector<Point>;
  using Triangle_range = std::vector<std::vector<std::size_t> >;

  {
    Grid grid { Point{-1., -1., -1.}, Point{1., 1., 1.}, std::array<std::size_t, 3>{2, 2, 2} };
    Values values { grid };
    Domain domain { grid, values };

    to_grid(grid, values, case_values);

    Point_range points;
    Triangle_range triangles;
    IS::marching_cubes<CGAL::Sequential_tag>(domain, iso, points, triangles, CGAL::parameters::use_topologically_correct_marching_cubes(true));

    assert(points.size() && triangles.size());
    assert(!has_duplicate_points(points, triangles));
    assert(!has_duplicate_polygons(points, triangles));
    assert(!has_isolated_vertices(points, triangles));

    assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles));
    Mesh m;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, m);

    assert(connected_components(m) == components);
    assert(euler_characteristic(m) == euler);
    assert(boundary_components(m) == boundaries);
  }

  {
    Grid grid_closed { Point{-1., -1., -1.}, Point{1., 1., 1.}, std::array<std::size_t, 3>{4, 4, 4} };
    Values values_closed { grid_closed };
    Domain domain_closed { grid_closed, values_closed };

    embed_in_positive_grid(grid_closed, values_closed, case_values);

    Point_range points;
    Triangle_range triangles;
    IS::marching_cubes<CGAL::Sequential_tag>(domain_closed, iso, points, triangles, CGAL::parameters::use_topologically_correct_marching_cubes(true));

    assert(points.size() && triangles.size());
    assert(!has_duplicate_points(points, triangles));
    assert(!has_duplicate_polygons(points, triangles));
    assert(!has_isolated_vertices(points, triangles));

    assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles));
    Mesh m;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, m);

    assert(!has_degenerate_faces(m));
    assert(is_manifold(m));
    assert(is_watertight(m));
  }
}

template <typename K>
void test_cube()
{
  using FT = typename K::FT;
  using Point = typename K::Point_3;
  using Vector = typename K::Vector_3;

  using Mesh = CGAL::Surface_mesh<Point>;
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Value_function_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;

  using Point_range = std::vector<Point>;
  using Polygon_range = std::vector<std::vector<std::size_t> >;

  using Mesh = CGAL::Surface_mesh<Point>;

  auto implicit_function = [](const Point& q) -> FT
  {
  // ---
    auto b1 = [](const Point& q) { return IS::Shapes::box<K>(Point(-1,-1,-1), Point(1,1,1), q); };
    auto s1 = [](const Point& q) { return IS::Shapes::sphere<K>(Point(0,0,0), 1, q); };
    auto cyl1 = [](const Point& q) { return IS::Shapes::infinite_cylinder<K>(Point(0,0,0), Vector(0,0,1), 0.5, q); };
    auto cyl2 = [](const Point& q) { return IS::Shapes::infinite_cylinder<K>(Point(0,0,0), Vector(0,1,0), 0.5, q); };
    auto cyl3 = [](const Point& q) { return IS::Shapes::infinite_cylinder<K>(Point(0,0,0), Vector(1,0,0), 0.5, q); };

    auto b1_is1 = [&](const Point& q) { return IS::Shapes::shape_intersection<K>(b1, s1, q); };
    auto cyl1_cyl2 = [&](const Point& q) { return IS::Shapes::shape_union<K>(cyl1, cyl2, q); };
    auto cyl1_cyl2_cyl3 = [&](const Point& q) { return IS::Shapes::shape_union<K>(cyl1_cyl2, cyl3, q); };

    auto b1_ms1_mcyl1_mcyl2 = [&](const Point& q) { return IS::Shapes::shape_difference<K>(b1_is1, cyl1_cyl2_cyl3, q); };

    return b1_ms1_mcyl1_mcyl2(q);
  };

  const CGAL::Bbox_3 bbox = {-2., -2., -2., 2., 2., 2.};
  const Vector spacing { 0.05, 0.05, 0.05 };
  const FT isovalue = 0;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running TMC (Implicit function) with isovalue = " << isovalue << std::endl;
  std::cout << "Kernel: " << typeid(K).name() << std::endl;

  Grid grid { bbox, spacing };
  Values values {implicit_function , grid };
  Domain domain { grid, values };

  Mesh debug_grid;
  auto debug_grid_creator = [&](const typename Domain::cell_descriptor& c)
    {
      std::vector<typename Mesh::Vertex_index> cell_vertices;
      for (const auto& v : domain.cell_vertices(c)) {
        cell_vertices.push_back(debug_grid.add_vertex(domain.point(v)));
      }
      debug_grid.add_face(cell_vertices[6], cell_vertices[2], cell_vertices[0], cell_vertices[4]);
      debug_grid.add_face(cell_vertices[1], cell_vertices[3], cell_vertices[7], cell_vertices[5]);
      debug_grid.add_face(cell_vertices[0], cell_vertices[1], cell_vertices[5], cell_vertices[4]);
      debug_grid.add_face(cell_vertices[6], cell_vertices[7], cell_vertices[3], cell_vertices[2]);
      debug_grid.add_face(cell_vertices[2], cell_vertices[3], cell_vertices[1], cell_vertices[0]);
      debug_grid.add_face(cell_vertices[4], cell_vertices[5], cell_vertices[7], cell_vertices[6]);
    };
    domain.template for_each_cell<CGAL::Sequential_tag>(debug_grid_creator);
  CGAL::IO::write_OFF("debug_grid.off", debug_grid);

  Point_range points;
  Polygon_range polygons;
  IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, isovalue, points, polygons, CGAL::parameters::use_topologically_correct_marching_cubes(true));

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #polygons: " << polygons.size() << std::endl;

  CGAL::IO::write_polygon_soup("test_cube.off", points, polygons, CGAL::parameters::stream_precision(17));
}

int main(int, char**)
{
  using K = CGAL::Simple_cartesian<double>;
  using FT = typename K::FT;

  std::array<FT, 8> case_values;
  FT iso = 0;

  // this is what the algorithm is built for, this works
  iso = generate_predefined_inner_ambiguity(AmbiguousCase::MC_13_TUNNEL, case_values);
  compare_tmc_mc_trilinear<K>(case_values, iso);  // export meshes to compare tmc and mc to the high resolution trilinear interpolation
  assert_tmc<K>(case_values, iso, 2, 1, 3);  // check if the mesh created by tmc has the expected topological properties

  // issue with degenerated hyperbolas
  iso = generate_predefined_singular(SingularCase::SIMPLE, case_values);
  compare_tmc_mc_trilinear<K>(case_values, iso);

  iso = generate_predefined_singular(SingularCase::TUNNEL, case_values);
  compare_tmc_mc_trilinear<K>(case_values, iso);

  // issue with edge with isovalue
  // the simple example still works
  iso = generate_predefined_iso_edge(IsoEdgeCase::SIMPLE, case_values);
  compare_tmc_mc_trilinear<K>(case_values, iso);
  assert_tmc<K>(case_values, iso, 1, 1, 1);

  // this does not properly work
  iso = generate_predefined_iso_edge(IsoEdgeCase::DOUBLE_SINGULAR, case_values);
  compare_tmc_mc_trilinear<K>(case_values, iso);
  // assert_tmc<K>(case_values, iso, 1, 1, 1);

  // this cetagory includes cases that are almost always non-manifold
  iso = generate_predefined_plane(PlaneCase::CROSS, case_values);
  // compare_tmc_mc_trilinear<K>(case_values, iso);

  // TODO: more examples still to test...

  // this is the example, where the issue was originally discovered
  // TODO: find which of the cases occur in this example and fail
  test_cube<K>();

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}

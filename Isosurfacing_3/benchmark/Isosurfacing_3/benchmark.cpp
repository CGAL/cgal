#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/internal/Isosurfacing_domain_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_gradients_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>

#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Image_3.h>
#include <CGAL/Isosurfacing_3/IO/Image_3.h>
#include <CGAL/Real_timer.h>

#include <cmath>
#include <iostream>
#include <limits>

#ifndef M_PI
# define M_PI 3.141592653589793238462643383279502884L
#endif

template <class GeomTraits>
struct Sphere_value
{
  using FT = typename GeomTraits::FT;
  using Point = typename GeomTraits::Point_3;

  FT operator()(const Point& point) const
  {
    return CGAL::approximate_sqrt((point - CGAL::ORIGIN).squared_length());
  }
};

template <class GeomTraits>
struct Sphere_gradient
{
  using FT = typename GeomTraits::FT;
  using Point = typename GeomTraits::Point_3;
  using Vector = typename GeomTraits::Vector_3;

  Vector operator()(const Point& point) const
  {
    Vector g = point - CGAL::ORIGIN;
    return g / CGAL::approximate_sqrt(g.squared_length());
  }
};

template <class GeomTraits>
struct Implicit_sphere
{
  using FT = typename GeomTraits::FT;
  using Vector = typename GeomTraits::Vector_3;

  using Grid = CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>;
  using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
  using Gradients = CGAL::Isosurfacing::Gradient_function_3<Grid>;
  using Domain = CGAL::Isosurfacing::internal::Isosurfacing_domain_3<Grid, Values, Gradients>;

  Implicit_sphere(const std::size_t N)
    : res(2. / N, 2. / N, 2. / N),
      grid { CGAL::Bbox_3 {-1, -1, -1, 1, 1, 1}, CGAL::make_array<std::size_t>(N, N, N) },
      values { Sphere_value<GeomTraits>{}, grid },
      gradients { Sphere_gradient<GeomTraits>{}, grid }
  { }

  Domain domain() const
  {
    return { grid, values, gradients };
  }

  FT iso() const
  {
    return 0.8;
  }

private:
  Vector res;
  Grid grid;
  Values values;
  Gradients gradients;
};

template <class GeomTraits>
struct IWPValue
{
  using FT = typename GeomTraits::FT;
  using Point = typename GeomTraits::Point_3;

  FT operator()(const Point& point) const
  {
    const FT alpha = 5.01;
    // const FT alpha = 1.01;

    const FT x = alpha * (point.x() + 1) * M_PI;
    const FT y = alpha * (point.y() + 1) * M_PI;
    const FT z = alpha * (point.z() + 1) * M_PI;

    return cos(x)*cos(y) + cos(y)*cos(z) + cos(z)*cos(x) - cos(x)*cos(y)*cos(z); // isovalue = 0
  }
};

template <class GeomTraits>
struct IWPGradient
{
  using FT = typename GeomTraits::FT;
  using Point = typename GeomTraits::Point_3;
  using Vector = typename GeomTraits::Vector_3;

  Vector operator()(const Point& point) const
  {
    const FT alpha = 5.01;
    // const FT alpha = 1.01;

    const FT x = alpha * (point.x() + 1) * M_PI;
    const FT y = alpha * (point.y() + 1) * M_PI;
    const FT z = alpha * (point.z() + 1) * M_PI;

    const FT gx = M_PI * alpha * sin(x) * (cos(y) * (cos(z) - 1.0) - cos(z));
    const FT gy = M_PI * alpha * sin(y) * (cos(x) * (cos(z) - 1.0) - cos(z));
    const FT gz = M_PI * alpha * sin(z) * (cos(x) * (cos(y) - 1.0) - cos(y));

    return { gx, gy, gz };
  }
};

template <class GeomTraits>
struct Implicit_iwp
{
  using FT = typename GeomTraits::FT;
  using Vector = typename GeomTraits::Vector_3;

  using Grid = CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>;
  using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
  using Gradients = CGAL::Isosurfacing::Gradient_function_3<Grid>;
  using Domain = CGAL::Isosurfacing::internal::Isosurfacing_domain_3<Grid, Values, Gradients>;

  Implicit_iwp(const std::size_t N)
    : res(2. / N, 2. / N, 2. / N)
  { }

  Domain domain() const
  {
    return { { {-1, -1, -1, 1, 1, 1}, res } , { IWPValue<GeomTraits>{} }, { IWPGradient<GeomTraits>{} } };
  }

  FT iso() const
  {
    return 0.;
  }

private:
  Vector res;
};

template <class GeomTraits>
struct Grid_sphere
{
  using Grid = CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>;
  using Values = CGAL::Isosurfacing::Interpolated_discrete_values_3<Grid>;
  using Gradients = CGAL::Isosurfacing::Interpolated_discrete_gradients_3<Grid>;
  using Domain = CGAL::Isosurfacing::internal::Isosurfacing_domain_3<Grid, Values, Gradients>;

  using FT = typename GeomTraits::FT;
  using Point = typename GeomTraits::Point_3;

  Grid_sphere(const std::size_t N)
  {
    const CGAL::Bbox_3 bbox{-1., -1., -1., 1., 1., 1.};
    grid = Grid { bbox, CGAL::make_array<std::size_t>(N, N, N) };

    values = { grid };
    gradients = { grid };

    const FT resolution = 2.0 / N;

    Sphere_value<GeomTraits> sphere_val;
    Sphere_gradient<GeomTraits> sphere_grad;

    for(std::size_t x = 0; x < grid.xdim(); x++)
    {
      const FT xp = x * resolution - 1.0;
      for(std::size_t y = 0; y < grid.ydim(); y++)
      {
        const FT yp = y * resolution - 1.0;
        for(std::size_t z = 0; z < grid.zdim(); z++)
        {
          const FT zp = z * resolution - 1.0;

          values(x, y, z) = sphere_val(Point(xp, yp, zp));
          gradients(x, y, z) = sphere_grad(Point(xp, yp, zp));
        }
      }
    }
  }

  Domain domain() const
  {
    return { grid, values, gradients };
  }

  typename GeomTraits::FT iso() const
  {
    return 0.8;
  }

private:
  Grid grid;
  Values values;
  Gradients gradients;
};

template <class GeomTraits>
struct Skull_image
{
  using Grid = CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>;
  using Values = CGAL::Isosurfacing::Interpolated_discrete_values_3<Grid>;
  using Gradients = CGAL::Isosurfacing::Finite_difference_gradient_3<GeomTraits>;
  using Domain = CGAL::Isosurfacing::internal::Isosurfacing_domain_3<Grid, Values, Gradients>;

  Skull_image(const std::size_t N)
  {
    const std::string fname = CGAL::data_file_path("images/skull_2.9.inr");
    CGAL::Image_3 image;
    if(!image.read(fname))
      std::cerr << "Error: Cannot read file " << fname << std::endl;

    Grid grid;
    Values values { grid };
    if(!CGAL::Isosurfacing::IO::read_Image_3(image, grid, values))
      std::cerr << "Error: Cannot convert image to Cartesian grid" << std::endl;

    gradients = { values };
  }

  Domain domain() const
  {
    return { grid, values, gradients };
  }

  typename GeomTraits::FT iso() const
  {
    return 2.9;
  }

private:
  Grid grid;
  Values values;
  Gradients gradients;
};

int main(int argc, char* argv[])
{
  std::size_t N = 100;

  const int argc_check = argc - 1;

  for(int i=1; i<argc; ++i)
  {
    if(!strcmp("-N", argv[i]) && i < argc_check)
      N = std::stoi(argv[++i]);
  }

#if defined KERNEL_SIMPLE_CARTESIAN_DOUBLE
  std::cout << "KERNEL_SIMPLE_CARTESIAN_DOUBLE" << std::endl;
  using Kernel = CGAL::Simple_cartesian<double>;
#elif defined KERNEL_SIMPLE_CARTESIAN_FLOAT
  std::cout << "KERNEL_SIMPLE_CARTESIAN_FLOAT" << std::endl;
  using Kernel = CGAL::Simple_cartesian<float>;
#elif defined KERNEL_CARTESIAN_DOUBLE
  std::cout << "KERNEL_CARTESIAN_DOUBLE" << std::endl;
  using Kernel = CGAL::Cartesian<double>;
#elif defined KERNEL_EPIC
  std::cout << "KERNEL_EPIC" << std::endl;
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
#else
  std::cout << "no kernel selected!" << std::endl;
  using Kernel = CGAL::Simple_cartesian<double>;
#endif

  using Point = Kernel::Point_3;

  using Point_range = std::vector<Point>;
  using Polygon_range = std::vector<std::vector<std::size_t> >;

#if defined SCENARIO_GRID_SPHERE
  std::cout << "SCENARIO_GRID_SPHERE" << std::endl;
  auto scenario = Grid_sphere<Kernel>(N);
#elif defined SCENARIO_IMPLICIT_SPHERE
  std::cout << "SCENARIO_IMPLICIT_SPHERE" << std::endl;
  auto scenario = Implicit_sphere<Kernel>(N);
#elif defined SCENARIO_IMPLICIT_IWP
  std::cout << "SCENARIO_IMPLICIT_IWP" << std::endl;
  auto scenario = Implicit_iwp<Kernel>(N);
#elif defined SCENARIO_SKULL_IMAGE
  std::cout << "SCENARIO_SKULL_IMAGE" << std::endl;
  auto scenario = Skull_image<Kernel>(N);
#else
  std::cout << "no scenario selected!" << std::endl;
  auto scenario = Implicit_sphere<Kernel>(N);
#endif

  Point_range points;
  Polygon_range polygons;

  CGAL::Real_timer timer;
  timer.start();

#if defined TAG_PARALLEL
  std::cout << "TAG_PARALLEL" << std::endl;
  using Tag = CGAL::Parallel_tag;
#elif defined TAG_SEQUENTIAL
  std::cout << "TAG_SEQUENTIAL" << std::endl;
  using Tag = CGAL::Sequential_tag;
#else
  std::cout << "no tag selected!" << std::endl;
  using Tag = CGAL::Sequential_tag;
#endif

#if defined ALGO_MARCHING_CUBES
  std::cout << "ALGO_MARCHING_CUBES" << std::endl;
  CGAL::Isosurfacing::marching_cubes<Tag>(scenario.domain(), scenario.iso(), points, polygons, false);
#elif defined ALGO_DUAL_CONTOURING
  std::cout << "ALGO_DUAL_CONTOURING" << std::endl;
  CGAL::Isosurfacing::dual_contouring<Tag>(scenario.domain(), scenario.iso(), points, polygons);
#else
  std::cout << "no algorithm selected!" << std::endl;
  CGAL::Isosurfacing::marching_cubes<Tag>(scenario.domain(), scenario.iso(), points, polygons);
#endif

  timer.stop();

  if(points.size() > std::numeric_limits<std::size_t>::max() - 2)
    std::cout << "This should never print and only prevents optimizations" << std::endl;

  std::cout << "internal timer: " << timer.time() << std::endl;
  std::cout << "internal polygons: " << polygons.size() << std::endl;
  std::cout << "internal points: " << points.size() << std::endl;
}

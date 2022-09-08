#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_grid_domain.h>
#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Implicit_domain.h>
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <math.h>

#include <iostream>
#include <limits>

#include "CLI11.hpp"
#include "Timer.h"

template <class GeomTraits>
struct SphereValue {
    typename GeomTraits::FT operator()(const typename GeomTraits::Point_3& point) const {
        return CGAL::approximate_sqrt((point - CGAL::ORIGIN).squared_length());
    }
};

template <class GeomTraits>
struct SphereGradient {
    typename GeomTraits::Vector_3 operator()(const typename GeomTraits::Point_3& point) const {
        const typename GeomTraits::Vector_3 g = point - CGAL::ORIGIN;
        return g / CGAL::approximate_sqrt(g.squared_length());
    }
};

template <class GeomTraits>
struct Implicit_sphere {

    typedef CGAL::Isosurfacing::Implicit_domain<GeomTraits, SphereValue<GeomTraits>, SphereGradient<GeomTraits>> Domain;
    typedef typename GeomTraits::Vector_3 Vector;

    Domain domain(const std::size_t N) const {

        const Vector res(2.0 / N, 2.0 / N, 2.0 / N);
        return Domain({-1, -1, -1, 1, 1, 1}, res, SphereValue<GeomTraits>(), SphereGradient<GeomTraits>());
    }

    typename GeomTraits::FT iso() const {
        return 0.8;
    }
};

template <class GeomTraits>
struct Grid_sphere {

    typedef CGAL::Isosurfacing::Cartesian_grid_domain<GeomTraits> Domain;
    typedef CGAL::Cartesian_grid_3<GeomTraits> Grid;
    typedef typename GeomTraits::FT FT;
    typedef typename GeomTraits::Point_3 Point;

    Domain domain(const std::size_t N) const {

        const FT resolution = 2.0 / N;
        SphereValue<GeomTraits> val;
        SphereGradient<GeomTraits> grad;

        Grid grid(N, N, N, {-1, -1, -1, 1, 1, 1});

        for (std::size_t x = 0; x < grid.xdim(); x++) {
            const FT xp = x * resolution - 1.0;

            for (std::size_t y = 0; y < grid.ydim(); y++) {
                const FT yp = y * resolution - 1.0;

                for (std::size_t z = 0; z < grid.zdim(); z++) {
                    const FT zp = z * resolution - 1.0;

                    grid.value(x, y, z) = val(Point(xp, yp, zp));
                    grid.gradient(x, y, z) = grad(Point(xp, yp, zp));
                }
            }
        }
        return Domain(grid);
    }

    typename GeomTraits::FT iso() const {
        return 0.8;
    }
};

template <class GeomTraits>
struct Skull_image {

    typedef CGAL::Isosurfacing::Cartesian_grid_domain<GeomTraits> Domain;
    typedef CGAL::Cartesian_grid_3<GeomTraits> Grid;

    Domain domain(const std::size_t N) const {

        const std::string fname = "../data/skull_2.9.inr";
        CGAL::Image_3 image;
        if (!image.read(fname)) {
            std::cerr << "Error: Cannot read file " << fname << std::endl;
            return EXIT_FAILURE;
        }

        const Grid grid(image);
        return Domain(grid);
    }

    typename GeomTraits::FT iso() const {
        return 2.9;
    }
};

int main(int argc, char* argv[]) {
    CLI::App app{"Isosurfacing benchmarks"};

    std::size_t N = 2;
    app.add_option("-N", N, "Grid size");

    CLI11_PARSE(app, argc, argv);

#if defined KERNEL_SIMPLE_CARTESIAN_DOUBLE
    std::cout << "KERNEL_SIMPLE_CARTESIAN_DOUBLE" << std::endl;
    typedef CGAL::Simple_cartesian<double> Kernel;
#elif defined KERNEL_SIMPLE_CARTESIAN_FLOAT
    std::cout << "KERNEL_SIMPLE_CARTESIAN_FLOAT" << std::endl;
    typedef CGAL::Simple_cartesian<float> Kernel;
#elif defined KERNEL_CARTESIAN_DOUBLE
    std::cout << "KERNEL_CARTESIAN_DOUBLE" << std::endl;
    typedef CGAL::Cartesian<double> Kernel;
#elif defined KERNEL_EPIC
    std::cout << "KERNEL_EPIC" << std::endl;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
#else
    std::cout << "no kernel selected!" << std::endl;
    typedef CGAL::Simple_cartesian<double> Kernel;
#endif

    typedef Kernel::Point_3 Point;

    typedef std::vector<Point> Point_range;
    typedef std::vector<std::vector<std::size_t>> Polygon_range;

#if defined SCENARIO_GRID_SPHERE
    std::cout << "SCENARIO_GRID_SPHERE" << std::endl;
    auto scenario = Grid_sphere<Kernel>();
#elif defined SCENARIO_IMPLICIT_SPHERE
    std::cout << "SCENARIO_IMPLICIT_SPHERE" << std::endl;
    auto scenario = Implicit_sphere<Kernel>();
#elif defined SCENARIO_SKULL_IMAGE
    std::cout << "SCENARIO_SKULL_IMAGE" << std::endl;
    auto scenario = Skull_image<Kernel>();
#else
    std::cout << "no scenario selected!" << std::endl;
    auto scenario = Implicit_sphere<Kernel>();
#endif

    Point_range points;
    Polygon_range polygons;

    ScopeTimer timer;

#if defined TAG_PARALLEL
    std::cout << "TAG_PARALLEL" << std::endl;
    typedef CGAL::Parallel_tag Tag;
#elif defined TAG_SEQUENTIAL
    std::cout << "TAG_SEQUENTIAL" << std::endl;
    typedef CGAL::Sequential_tag Tag;
#else
    std::cout << "no tag selected!" << std::endl;
    typedef CGAL::Sequential_tag Tag;
#endif

#if defined ALGO_MARCHING_CUBES
    std::cout << "ALGO_MARCHING_CUBES" << std::endl;
    CGAL::Isosurfacing::make_triangle_mesh_using_marching_cubes<Tag>(scenario.domain(N), scenario.iso(), points,
                                                                     polygons);
#elif defined ALGO_DUAL_CONTOURING
    std::cout << "ALGO_DUAL_CONTOURING" << std::endl;
    CGAL::Isosurfacing::make_quad_mesh_using_dual_contouring<Tag>(scenario.domain(N), scenario.iso(), points, polygons);
#else
    std::cout << "no algorithm selected!" << std::endl;
    CGAL::Isosurfacing::make_triangle_mesh_using_marching_cubes<Tag>(scenario.domain(N), scenario.iso(), points,
                                                                     polygons);
#endif

    const int64_t ms = timer.stop();

    if (points.size() > std::numeric_limits<std::size_t>::max() - 2) {
        std::cout << "This should never print and only prevents optimizations" << std::endl;
    }

    std::cout << "internal timer: " << ms << std::endl;
}

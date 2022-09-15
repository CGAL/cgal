#ifndef CGAL_ISOSURFACING_TEST_UTIL_H
#define CGAL_ISOSURFACING_TEST_UTIL_H

#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_grid_domain.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Implicit_domain.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Surface_mesh.h>
#include <tbb/concurrent_vector.h>

#include <cassert>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Point_3 Point;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

namespace PMP = CGAL::Polygon_mesh_processing;

bool has_duplicate_points(Point_range points, Polygon_range polygons) {
    return PMP::merge_duplicate_points_in_polygon_soup(points, polygons) != 0;
}

bool has_duplicate_polygons(Point_range points, Polygon_range polygons) {
    return PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons) != 0;
}

bool has_isolated_vertices(Point_range points, Polygon_range polygons) {
    return PMP::remove_isolated_points_in_polygon_soup(points, polygons) != 0;
}

bool is_polygon_mesh(const Polygon_range& polygons) {
    return PMP::is_polygon_soup_a_polygon_mesh(polygons);
}

Mesh to_mesh(const Point_range& points, const Polygon_range& polygons) {
    Mesh m;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, m);
    return m;
}

bool is_manifold(Mesh& m) {
    return PMP::duplicate_non_manifold_vertices(m, CGAL::parameters::dry_run(true)) == 0;
}

bool has_degenerate_faces(Mesh& m) {
    return PMP::remove_connected_components_of_negligible_size(
               m, CGAL::parameters::dry_run(true).area_threshold(std::numeric_limits<FT>::epsilon())) != 0;
}

bool check_mesh_distance(const Mesh& m0, const Mesh& m1) {
    auto dist = PMP::approximate_Hausdorff_distance<CGAL::Sequential_tag>(
        m0, m1, CGAL::parameters::number_of_points_per_area_unit(4000));
    std::cout << dist << std::endl;
    return true;
}

#endif  // CGAL_ISOSURFACING_TEST_UTIL_H
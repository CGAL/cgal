// This code was first submitted in an issue:
//    https://github.com/CGAL/cgal/issues/4025
// and then rewrote a lot, keeping the observed behavior.

#include <CGAL/Base_with_time_stamp.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <string_view>
#include <iostream>
#include <iterator>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polygon_2 = CGAL::Polygon_2<K>;
using Itag_ = CGAL::Exact_intersections_tag;
using Vb = CGAL::Base_with_time_stamp<CGAL::Triangulation_vertex_base_2<K>>;
using Cb = CGAL::Base_with_time_stamp<CGAL::Constrained_triangulation_face_base_2<K>>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Cb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag_>;
using CDTP = CGAL::Constrained_triangulation_plus_2<CDT>;

using Point = CDTP::Point;
using Vertex_handle = CDTP::Vertex_handle;
using Constraint_id = CDTP::Constraint_id;
using Vertices_in_constraint_iterator = CDTP::Vertices_in_constraint_iterator;

auto nb_of_vertices(CDTP &cdtp, Constraint_id id)
{
    return static_cast<std::size_t>(std::distance(cdtp.vertices_in_constraint_begin(id),
                                                  cdtp.vertices_in_constraint_end(id)));
}

auto value_check_expected =
    [](auto&& value, [[maybe_unused]] const auto& expected) -> decltype(auto)
    {
        assert(value == expected);
        return std::forward<decltype(value)>(value);
    };

auto oformat =
    [](Vertex_handle vh)
    {
        return CGAL::IO::oformat(vh, CGAL::With_point_tag{});
    };

int main()
{
    CDTP cdtp;

    auto print_cdtp = [&cdtp](std::string_view msg)
    {
        std::cout << msg << std::endl;
        cdtp.print_hierarchy(std::cout);
    };

    const std::array<Point, 6> collinear_points = {
        Point(0,0), Point(0,1), Point(0,2), Point(0,3), Point(0,4), Point(0,5)
    };

    const std::array<Point, 6> non_collinear_points = {
        Point(1,0), Point(2,1), Point(4,2), Point(2,3), Point(4,4), Point(1,5)
    };

    Constraint_id collinear_cid = cdtp.insert_constraint(collinear_points.begin(), collinear_points.end());
    Constraint_id non_collinear_cid = cdtp.insert_constraint(non_collinear_points.begin(), non_collinear_points.end());

    print_cdtp("Initial state");

    Vertices_in_constraint_iterator vertex_it = std::next(cdtp.vertices_in_constraint_begin(collinear_cid), 2);
    [[maybe_unused]] auto next_it = std::next(vertex_it);
    std::cout << "\n-> attempt to remove vertex " << oformat(*vertex_it) << std::endl;
    vertex_it = cdtp.remove_vertex_from_constraint(collinear_cid, vertex_it);
    std::cout << "   cdtp.remove_vertex_from_constraint(collinear_cid, vertex_it) returned the vertex "
              << oformat(*vertex_it) << std::endl;
    assert(vertex_it == next_it);

    print_cdtp("\nAfter removing third vertex from the collinear constraint");

    // The first constraint (ID `collinear_cid`) is collinear. `cdtp.remove_vertex_from_constraint`
    // cannot remove the third vertex from it, because it is collinear with the triangulation vertex
    // with the point (0, 2).

    std::cout << "\nnumber of subconstraints: "
              << value_check_expected(cdtp.number_of_subconstraints(), 10U) << std::endl;
    std::cout << "number of constraints: "
              << value_check_expected(cdtp.number_of_constraints(), 2U) << std::endl;
    std::cout << "number of vertex in collinear constraint: "
              << value_check_expected(nb_of_vertices(cdtp, collinear_cid), 6U) << std::endl;

    vertex_it = std::next(cdtp.vertices_in_constraint_begin(non_collinear_cid), 2);
    next_it = std::next(vertex_it);
    std::cout << "\n-> attempt to remove vertex " << oformat(*vertex_it) << std::endl;
    vertex_it = cdtp.remove_vertex_from_constraint(non_collinear_cid, vertex_it);
    std::cout << "   cdtp.remove_vertex_from_constraint(non_collinear_cid, vertex_it) returned the vertex "
              << oformat(*vertex_it) << std::endl;
    assert(vertex_it == next_it);

    print_cdtp("\nAfter removing third vertex from the non-collinear constraint");

    std::cout << "number of subconstraints: "
              << value_check_expected(cdtp.number_of_subconstraints(), 9U) << std::endl;
    std::cout << "number of constraints: "
              << value_check_expected(cdtp.number_of_constraints(), 2U) << std::endl;
    std::cout << "number of vertex in collinear constraint: "
              << value_check_expected(nb_of_vertices(cdtp, non_collinear_cid), 5U) << std::endl;

    return 0;
}

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

template <typename V, typename E>
decltype(auto) value_check_expected(V&& value, [[maybe_unused]] const E& expected)
{
  assert(value == expected);
  return std::forward<V>(value);
};

auto oformat(Vertex_handle vh)
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


    // NOW test another scenario, that has nothing to do with the issue #4025
    cdtp.clear();
    print_cdtp("\nAfter clearing the constrained triangulation");


    // Let's insert a constraint with a loop
    //                  (1,1)
    //                  /  |
    //                 /   |
    //  start-->(0,0) X-->(1,0)
    //               /
    //              /
    //           (-1,-1)
    const std::array<Point, 4> looping_cid = {
        Point(0,0), Point(1,0), Point(1,1), Point(-1,-1)
    };

    cdtp.insert_constraint(looping_cid.begin(), looping_cid.end());
    print_cdtp("\nAfter inserting a looping constraint");
    std::cout << "\nnumber of subconstraints: "
              << value_check_expected(cdtp.number_of_subconstraints(), 4U) << std::endl;

    // NOW test another scenario
    cdtp.clear();
    print_cdtp("\nAfter clearing the constrained triangulation");
    // Let's insert a constraint with identical sub-constraints
    //               (1,1)
    //              /   |
    //             /    |
    // start-->(0,0)-->(1,0)--->(3,0)
    const std::array<Point, 5> overlaping_cid = {Point(0, 0), Point(1, 0), Point(1, 1),
                                                 Point(0, 0), Point(3, 0)};
    cdtp.insert_constraint(overlaping_cid.begin(), overlaping_cid.end());
    print_cdtp("\nAfter inserting a constraint with overlapping subconstraints");
    std::cout << "\nnumber of subconstraints: "
              << value_check_expected(cdtp.number_of_subconstraints(), 4U)
              << "\ncdtp.subconstraints.size(): "
              << value_check_expected(cdtp.subconstraints().size(), 4U) << std::endl;

    // NOW test another scenario
    cdtp.clear();
    print_cdtp("\nAfter clearing the constrained triangulation");
    // Let's insert two constraints with four points each and one common segment in the middle
    //     start-->(0,1)             (3,1)
    //                 \             /
    //                  \           /
    // start-->(0,0)--->(1,0)===>(2,0)--->(3,0)
    const std::array<Point, 4> first_cid =  {Point(0, 0), Point(1, 0), Point(2, 0), Point(3, 0)};
    const std::array<Point, 4> second_cid = {Point(0, 1), Point(1, 0), Point(2, 0), Point(3, 1)};
    cdtp.insert_constraint(first_cid.begin(), first_cid.end());
    cdtp.insert_constraint(second_cid.begin(), second_cid.end());
    print_cdtp("\nAfter inserting two constraints with a common segment in the middle");
    std::cout << "\nnumber of subconstraints: "
              << value_check_expected(cdtp.number_of_subconstraints(), 5U)
              << "\ncdtp.subconstraints.size(): "
              << value_check_expected(cdtp.subconstraints().size(), 5U) << std::endl;

    return 0;
}

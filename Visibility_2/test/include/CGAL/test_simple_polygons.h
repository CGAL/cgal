/*
 * Author: Francisc Bungiu 
 * E-mail: fbungiu@gmail.com
 */

#ifndef CGAL_TEST_SIMPLE_POLYGONS_H
#define CGAL_TEST_SIMPLE_POLYGONS_H

namespace CGAL {

template < class _Visibility_2, class _Arrangement_2 >
bool simple_polygon_test_case_1() {

	typedef _Visibility_2									Visibility_2;
	typedef _Arrangement_2 									Arrangement_2;
	typedef typename Arrangement_2::Geometry_traits_2		Geometry_traits_2;
	typedef typename Geometry_traits_2::Point_2				Point_2;

	Visibility_2 visibility;

	// First read arrangement 
	Arrangement_2 arr, out_arr, correct_out_arr;
    std::ifstream input("./data/simple_polygon_test_case_1.in");
    CGAL::create_arrangement_from_input<Arrangement_2>(arr, input);

	// Read query point from file
	double x, y;
    input >> x >> y;
    Point_2 query_pt(x, y);

    std::ifstream correct_output("./data/simple_polygon_test_case_1.out");
    CGAL::create_arrangement_from_input<Arrangement_2>(correct_out_arr, correct_output);

	typename Arrangement_2::Face_const_iterator fit;

    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (!fit->is_unbounded()) {
            visibility.visibility_region(query_pt, fit, out_arr);
        }
    }

    return CGAL::test_are_equal<Arrangement_2>(correct_out_arr, out_arr);
}

} // end namespace CGAL

#endif
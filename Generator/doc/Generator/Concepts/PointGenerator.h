/*!
\ingroup PkgGeneratorsConcepts
\cgalConcept

The concept `PointGenerator` defines the requirements for a point generator, 
which can be used in places where input iterators are called for. 

\cgalHasModel `CGAL::Random_points_in_ball_d<Point_d,>`
\cgalHasModel `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\cgalHasModel `CGAL::Random_points_in_square_2<Point_2, Creator>`
\cgalHasModel `CGAL::Random_points_in_triangle_2<Point_2, Creator>`
\cgalHasModel `CGAL::Random_points_on_circle_2<Point_2, Creator>`
\cgalHasModel `CGAL::Random_points_on_segment_2<Point_2, Creator>`
\cgalHasModel `CGAL::Random_points_on_square_2<Point_2, Creator>`
\cgalHasModel `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\cgalHasModel `CGAL::Random_points_in_cube_d<Point_d>`
\cgalHasModel `CGAL::Random_points_in_sphere_3<Point_3, Creator>`
\cgalHasModel `CGAL::Random_points_in_triangle_3<Point_3, Creator>`
\cgalHasModel `CGAL::Random_points_in_tetrahedron_3<Point_3, Creator>`
\cgalHasModel `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\cgalHasModel `CGAL::Random_points_on_sphere_d<Point_d>`

*/

class PointGenerator {
public:

/// \name Types 
/// @{

/*!
the type of point being generated. 
*/ 
typedef unspecified_type value_type; 

/// @} 

/// \name Operations 
/// @{

/*!
returns an absolute bound for 
the coordinates of all generated points. 
*/ 
double range() const; 

/// @}

}; /* end PointGenerator */

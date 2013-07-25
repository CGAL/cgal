
namespace CGAL {

/*!
\ingroup PkgInterpolation2SurfaceNeighbor

`Voronoi_intersection_2_traits_3` is a model for the concept 
`RegularTriangulationTraits_2` and 
`InterpolationTraits`. It can be used to instantiate the 
geometric traits class of a two-dimensional regular triangulation. 
A three-dimensional plane is defined by a point and a vector that 
are members of the traits class. The triangulation is defined on `3D` 
points. It is the regular triangulation of the input points 
projected onto the plane and each weighted with the negative squared 
distance of the input point to the plane. It can be shown that it is 
dual to the power diagram obtained by intersecting the 
three-dimensional Voronoi diagram of the input points with the 
plane. All predicates and constructions used in the computation of 
the regular triangulation are formulated on the three dimensional 
points without explicitly constructing the projected points and the 
weights. This reduces the arithmetic demands. The traits class is 
templated by a kernel class `K`. 

\cgalModels `RegularTriangulationTraits_2`

\sa `RegularTriangulationTraits_2` 
\sa `CGAL::Regular_triangulation_2<Gt, Tds>` 
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa PkgInterpolationSurfaceNeighborCoordinates3

*/
template< typename K >
class Voronoi_intersection_2_traits_3 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef K::RT Weight; 

/*!

*/ 
typedef K::FT FT; 

/*!

*/ 
typedef K::Point_3 Point_2; 

/*!

*/ 
typedef K::Segment_3 Segment_2; 

/*!

*/ 
typedef K::Triangle_3 Triangle_2; 

/*!

*/ 
typedef K::Line_3 Line_2; 

/*!

*/ 
typedef K::Ray_3 Ray_2; 

/*!

*/ 
typedef K::Vector_3 Vector_2; 

/*!

*/ 
typedef K::Construct_triangle_3 Construct_triangle_2; 

/*!

*/ 
typedef K::Construct_ray_3 Construct_ray_2; 

/*!

*/ 
typedef K::Compare_distance_3 Compare_distance_2; 

/*!
An instance of this function object class computes the square 
root of the result of `K::Compute_squared_area_3`. 
If the number type `FT` does not support the square root 
operation, the result is cast to `double` 
before computing the square root. 
*/ 
Compute_area_3<Rep> Compute_area_2; 

/*!

*/ 
typedef Orientation_with_normal_plane_2_3<Rep> Orientation_2; 

/*!

*/ 
typedef Side_of_plane_centered_sphere_2_3<Point_2> Power_test_2; 

/*!

*/ 
typedef Construct_plane_centered_circumcenter_3<Point_2> 
Construct_weighted_circumcenter_2; 

/*!

*/ 
typedef Compare_first_projection_3<Point_2> Compare_x_2; 

/*!

*/ 
typedef Compare_second_projection_3<Point_2> Compare_y_2; 

/// @} 

/// \name Creation 
/// @{

/*!
The plane associated to the traits class contains `point` and has as normal vector `normal`. 
*/ 
Voronoi_intersection_2_traits_3(const typename K::Point_3& point=typename K::Point_3(), const typename K::Vector_3& normal=NULL_VECTOR); 

/// @}

}; /* end Voronoi_intersection_2_traits_3 */
} /* end namespace CGAL */

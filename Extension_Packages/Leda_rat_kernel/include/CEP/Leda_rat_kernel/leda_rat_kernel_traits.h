#ifndef CEP_LEDA_RAT_KERNEL_TRAITS_H
#define CEP_LEDA_RAT_KERNEL_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>


// I had to add the log function for LEDA integers to CGAL namespace 
// for the VC compiler (otherwise log in CGAL/leda_integer.h causes trouble ...)
#if defined(LEDA_NAMESPACE) && defined(_MSC_VER)
#include <CGAL/LEDA_basic.h>
#include <LEDA/integer.h>

namespace CGAL {
int log(const leda_integer& i)
{ return leda::log(i); }
}
#endif


// this is our "help-kernel" used
// for computations and intersections
// this garanties same results of the LEDA and CGAL kernels ...
#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <LEDA/basic.h>

// if we use a LEDA version without namespaces
// we have to define a few macros
#if !defined(LEDA_NAMESPACE)
#define LEDA_BEGIN_NAMESPACE
#define LEDA_END_NAMESPACE
#define LEDA_NAMESPACE_NAME
#endif

// 2d geometry ...
#include <LEDA/rat_point.h>
#include <LEDA/rat_segment.h>
#include <LEDA/rat_line.h>
#include <LEDA/rat_circle.h>
#include <LEDA/rat_ray.h>
#include <LEDA/rat_triangle.h>
#include <LEDA/rat_rectangle.h>
#include <LEDA/rat_vector.h>

// event support for CGAL 
#if defined(CGAL_GEOMETRY_EVENTS)
#include <CEP/Leda_rat_kernel/event.h>
#endif

// I had to write an own type for LEDA directions;
// this was necessary because in the Equal_2 kernel
// concept I had to provide operators taking a vector
// and direction in one class, so that "faking" the direction
// by using a LEDA rat_vector was not possible any longer ...
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/rat_direction.h>

#define CGAL_COMPATIBLE_CIRCLES
#define CGAL_COMPATIBLE_SPHERES

// do we want CGAL compatible circles ???
#if defined(CGAL_COMPATIBLE_CIRCLES)
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/cgal_rat_circle.h>
#endif

// do we want CGAL compatible spheres ???
#if defined(CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3)

#if defined(CGAL_COMPATIBLE_SPHERES)
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/cgal_d3_rat_sphere.h>
#endif

// 3d geometry ...
#include <LEDA/d3_rat_point.h>
#include <LEDA/d3_rat_plane.h>
#include <LEDA/d3_rat_segment.h>
#include <LEDA/d3_rat_triangle.h>
#include <LEDA/d3_rat_simplex.h>
#include <LEDA/d3_rat_sphere.h>
#include <LEDA/d3_rat_line.h>
#include <LEDA/d3_rat_ray.h>

#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/d3_rat_iso_cuboid.h>
#endif

#if defined(CGAL_PROVIDE_LEDA_PARTITION_TRAITS)
#include <CGAL/Polygon_2.h>
#include <list>
#endif

// -------------------------------------------------
//  ..... 2d kernel .....
// -------------------------------------------------
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/cgal_leda_conversion_2.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/constructions_2.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/intersections_2.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/computations_2.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/generalized_predicates_2.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/extension_predicates_2.h>

// -------------------------------------------------
//  ..... 3d kernel .....
// -------------------------------------------------
#if defined(CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3)
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/cgal_leda_conversion_3.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/constructions_3.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/intersections_3.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/computations_3.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/generalized_predicates_3.h>
#endif


CGAL_BEGIN_NAMESPACE

// kernel tag for the rational LEDA kernel traits

struct Leda_rational_kernel_tag {  };


class leda_rat_R {
public:
  typedef leda_rational             FT;
  typedef leda_integer              RT;
};

class leda_rat_kernel_traits {
public:

  // this is our "help kernel"
  // we use it for computing intersections
  // and for doing (a few) computations
  typedef CGAL::Homogeneous<leda_integer> HELP_KERNEL;
  
  typedef leda_rat_kernel_traits    Self;

  typedef leda_rat_R                R;
  
  //tag type ...
  typedef Leda_rational_kernel_tag  Kernel_tag;

  // types ...
  typedef leda_rational             FT;
  typedef leda_integer              RT;
  typedef leda_rat_point            Point_2;
  typedef leda_rat_vector           Vector_2;
  typedef LEDA_NAMESPACE_NAME::rat_direction           Direction_2;
  typedef leda_rat_line             Line_2;
  typedef leda_rat_ray              Ray_2;
  typedef leda_rat_segment          Segment_2;
  typedef leda_rat_triangle         Triangle_2;
  typedef leda_rat_rectangle        Iso_rectangle_2; 
  
#if defined(CGAL_COMPATIBLE_CIRCLES)
  typedef LEDA_NAMESPACE_NAME::cgal_rat_circle  Circle_2;
#else  
  // use LEDA rat_circles
  typedef leda_rat_circle           Circle_2;
#endif

  typedef CGAL::Object              Object_2;

  // typedef without _2 (for compatibility to old concepts)
  typedef leda_rational             Coord_type; // for alpha shapes ...
  typedef leda_rat_point            Point;
  typedef leda_rat_vector           Vector;
  typedef LEDA_NAMESPACE_NAME::rat_direction       Direction;
  typedef leda_rat_line             Line;
  typedef leda_rat_ray              Ray;
  typedef leda_rat_segment          Segment;
  typedef leda_rat_triangle         Triangle;
  typedef leda_rat_rectangle        Iso_rectangle;
#if defined(CGAL_COMPATIBLE_CIRCLES)
  typedef LEDA_NAMESPACE_NAME::cgal_rat_circle  Circle;
#else  
  // use LEDA rat_circles
  typedef leda_rat_circle           Circle;
#endif


// -----------------------------------------------------------
// support for 3d kernel traits ...
// -----------------------------------------------------------   

#if defined(CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3)
  typedef leda_d3_rat_point         Point_3;
  typedef leda_rat_vector           Vector_3;
  typedef LEDA_NAMESPACE_NAME::rat_direction     Direction_3;
  typedef LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid    Iso_cuboid_3;

  typedef leda_d3_rat_line          Line_3;
  typedef leda_d3_rat_ray           Ray_3;
  typedef leda_d3_rat_segment       Segment_3;
#if defined(CGAL_COMPATIBLE_SPHERES)
  typedef LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere  Sphere_3;
#else 
  typedef leda_d3_rat_sphere        Sphere_3;
#endif  
  typedef leda_d3_rat_plane         Plane_3;
  typedef leda_d3_rat_triangle      Triangle_3;
  typedef LEDA_NAMESPACE_NAME::d3_rat_simplex       Tetrahedron_3;
  typedef CGAL::Object              Object_3;
#endif  
// -----------------------------------------------------------   
 
#if defined(CGAL_PROVIDE_LEDA_PARTITION_TRAITS)
  typedef CGAL::Polygon_2<CGAL::leda_rat_kernel_traits, std::list<leda_rat_point> > Polygon_2;
#endif  
  
  // 2d construction types ...
  typedef Construct_leda_rat_point<Self>                          Construct_point_2;
  typedef Construct_leda_rat_vector<Self>                         Construct_vector_2;
  typedef Construct_leda_rat_direction<Self>                      Construct_direction_2;
  typedef Construct_leda_rat_direction<Self>                      Construct_direction_of_line_2;
  
  typedef Construct_leda_rat_segment<Self>                        Construct_segment_2;
  typedef Construct_leda_rat_line<Self>                           Construct_line_2;
  typedef Construct_leda_rat_ray<Self>                            Construct_ray_2;
  typedef Construct_leda_rat_circle<Self>                         Construct_circle_2;
  typedef Construct_leda_rat_triangle<Self>                       Construct_triangle_2;
  typedef Construct_leda_rat_rectangle<Self>                      Construct_iso_rectangle_2;
  typedef Construct_leda_rat_object                               Construct_object_2;
  typedef Construct_leda_rat_scaled_vector_2<Self>                Construct_scaled_vector_2;
  typedef Construct_leda_rat_translated_point_2<Self>             Construct_translated_point_2;
  typedef Construct_leda_rat_point_on_2<Self>                     Construct_point_on_2;
  typedef Construct_leda_rat_projected_point_2<Self>              Construct_projected_point_2;
  typedef Construct_leda_rat_vertex_2<Self>                       Construct_vertex_2;
  typedef Construct_leda_rat_supporting_line_2<Self>              Construct_supporting_line_2;
  typedef Construct_leda_rat_perpendicular_vector_2<Self>         Construct_perpendicular_vector_2;
  typedef Construct_leda_rat_perpendicular_direction_2<Self>      Construct_perpendicular_direction_2;
  typedef Construct_leda_rat_perpendicular_line_2<Self>           Construct_perpendicular_line_2;
  typedef Construct_leda_rat_midpoint_2<Self>                     Construct_midpoint_2;      
  typedef Construct_leda_rat_center_2<Self>                       Construct_center_2; 
  typedef Construct_leda_rat_centroid_2<Self>                     Construct_centroid_2;    
  typedef Construct_leda_rat_circumcenter_2<Self>                 Construct_circumcenter_2;
  typedef Construct_leda_rat_bisector_2<Self>                     Construct_bisector_2;
  typedef Construct_leda_rat_opposite_direction_2<Self>           Construct_opposite_direction_2;
  typedef Construct_leda_rat_opposite_segment_2<Self>             Construct_opposite_segment_2;      
  typedef Construct_leda_rat_opposite_ray_2<Self>                 Construct_opposite_ray_2;    
  typedef Construct_leda_rat_opposite_line_2<Self>                Construct_opposite_line_2;
  typedef Construct_leda_rat_opposite_triangle_2<Self>            Construct_opposite_triangle_2;
  typedef Construct_leda_rat_opposite_circle_2<Self>              Construct_opposite_circle_2;  
  typedef Construct_leda_rat_opposite_vector_2<Self>              Construct_opposite_vector_2;
  
  // 2d intersections and related things
  typedef Assign_leda_rat_2                                       Assign_2;
  typedef CGAL_intersect_leda_rat_2<Self,HELP_KERNEL>             Intersect_2;
    
  // 2d computations
  typedef CGAL_compute_leda_rat_squared_distance_2<HELP_KERNEL>  Compute_squared_distance_2;
  typedef CGAL_compute_leda_rat_y_at_x_2<HELP_KERNEL>            Compute_y_at_x_2;  
  typedef Compute_leda_rat_squared_length_2<Self>                Compute_squared_length_2;
  typedef Compute_leda_rat_squared_radius_2<Self>                Compute_squared_radius_2;
  typedef Compute_leda_rat_area_2<Self>                          Compute_area_2;    
    
  // 2d generalized predicates
  typedef Predicate_leda_rat_angle_2<Self>                        Angle_2;
  typedef Predicate_leda_rat_equal_2<Self>                        Equal_2;
  typedef Predicate_leda_rat_equal_x_2<Self>                      Equal_x_2;
  typedef Predicate_leda_rat_equal_y_2<Self>                      Equal_y_2;
  typedef Predicate_leda_rat_less_x_2<Self>                       Less_x_2;
  typedef Predicate_leda_rat_less_x_2<Self>                       Less_y_2;
  typedef Predicate_leda_rat_less_xy_2<Self>                      Less_xy_2;
  typedef Predicate_leda_rat_less_yx_2<Self>                      Less_yx_2;
  typedef Predicate_leda_rat_compare_x_2<Self>                    Compare_x_2;
  typedef Predicate_leda_rat_compare_x_at_y_2<Self>               Compare_x_at_y_2;
  typedef Predicate_leda_rat_compare_y_2<Self>                    Compare_y_2; 
  typedef Predicate_leda_rat_compare_xy_2<Self>                   Compare_xy_2;
  typedef Predicate_leda_rat_compare_y_at_x_2<Self>               Compare_y_at_x_2;
  typedef Predicate_leda_rat_compare_distance_2<Self>             Compare_distance_2;
  typedef Predicate_leda_rat_compare_angle_with_x_axis_2<Self>    Compare_angle_with_x_axis_2;
  
  // ----------------------------------------------------------------------------------------
  // Equal_xy_2 is undocumented, but in the interface macros ... 
  typedef Predicate_leda_rat_equal_xy_2<Self>                     Equal_xy_2;
  // ----------------------------------------------------------------------------------------
  
  typedef Predicate_leda_rat_compare_slope_2<Self>                Compare_slope_2;
  typedef Predicate_leda_rat_less_distance_to_point_2<Self>       Less_distance_to_point_2;  
  typedef Predicate_leda_rat_less_signed_distance_to_line_2<Self> Less_signed_distance_to_line_2;
  typedef Predicate_leda_rat_less_rotate_ccw_2<Self>              Less_rotate_ccw_2;
  typedef Predicate_leda_rat_leftturn_2<Self>                     Leftturn_2; 
  typedef Predicate_leda_rat_leftturn_2<Self>                     Left_turn_2;  
  typedef Predicate_leda_rat_collinear_2<Self>                    Collinear_2;
  typedef Predicate_leda_rat_orientation_2<Self>                  Orientation_2;
  typedef Predicate_leda_rat_side_of_oriented_circle_2<Self>      Side_of_oriented_circle_2;
  typedef Predicate_leda_rat_side_of_bounded_circle_2<Self>       Side_of_bounded_circle_2;
  typedef Predicate_leda_rat_is_horizontal_2<Self>                Is_horizontal_2;
  
  // additional  functor ...
  typedef Predicate_leda_rat_is_in_x_range_2<Self>                Is_in_x_range_2;
  
  typedef Predicate_leda_rat_is_vertical_2<Self>                  Is_vertical_2;
  typedef Predicate_leda_rat_is_degenerate_2<Self>                Is_degenerate_2;
  typedef Predicate_leda_rat_has_on_2<Self>                       Has_on_2;
  typedef Predicate_leda_rat_collinear_has_on_2<Self>             Collinear_has_on_2;
  typedef Predicate_leda_rat_has_on_bounded_side_2<Self>          Has_on_bounded_side_2;
  typedef Predicate_leda_rat_has_on_unbounded_side_2<Self>        Has_on_unbounded_side_2;
  typedef Predicate_leda_rat_has_on_boundary_2<Self>              Has_on_boundary_2;
  typedef Predicate_leda_rat_has_on_positive_side_2<Self>         Has_on_positive_side_2;
  typedef Predicate_leda_rat_has_on_negative_side_2<Self>         Has_on_negative_side_2;
  typedef Predicate_leda_rat_oriented_side_2<Self>                Oriented_side_2;
  typedef Predicate_leda_rat_bounded_side_2<Self>                 Bounded_side_2;
  typedef Predicate_leda_rat_are_ordered_along_line_2<Self>       Are_ordered_along_line_2;
  typedef Predicate_leda_rat_are_strictly_ordered_along_line_2<Self>   Are_strictly_ordered_along_line_2;
  typedef Predicate_leda_rat_collinear_are_ordered_along_line_2<Self>  Collinear_are_ordered_along_line_2;
  typedef Predicate_leda_rat_collinear_are_strictly_ordered_along_line_2<Self>  Collinear_are_strictly_ordered_along_line_2;
  typedef Predicate_leda_rat_counterclockwise_in_between_2<Self>  Counterclockwise_in_between_2;
  typedef Predicate_leda_rat_do_intersect_2<Self>                 Do_intersect_2;    

  // a number of extensions ...
  typedef Predicate_leda_rat_do_intersect_to_right_2<Self>        Do_intersect_to_right_2;
  typedef Predicate_leda_rat_do_intersect_to_left_2<Self>         Do_intersect_to_left_2;

// -----------------------------------------------------------
// support for 3d kernel traits ...
#if defined(CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3)
   // 3d constructions ...
   typedef Construct_leda_d3_rat_point<Self>                      Construct_point_3;
   typedef Construct_leda_d3_rat_vector<Self>                     Construct_vector_3;
   typedef Construct_leda_d3_rat_direction<Self>                  Construct_direction_3;
   
   // for 3d triangulations:
   typedef Construct_leda_d3_rat_direction<Self>                  Construct_direction_of_line_3;
   
   typedef Construct_leda_d3_rat_plane<Self>                      Construct_plane_3;
   typedef Construct_leda_d3_rat_iso_cuboid<Self>                 Construct_iso_cuboid_3;
   typedef Construct_leda_d3_rat_line<Self>                       Construct_line_3;
   typedef Construct_leda_d3_rat_ray<Self>                        Construct_ray_3;
   typedef Construct_leda_d3_rat_sphere<Self>                     Construct_sphere_3;     
   typedef Construct_leda_d3_rat_segment<Self>                    Construct_segment_3;
   typedef Construct_leda_d3_rat_simplex<Self>                    Construct_tetrahedron_3;
   typedef Construct_leda_d3_rat_triangle<Self>                   Construct_triangle_3;
   typedef Construct_leda_d3_rat_object                           Construct_object_3;
   typedef Construct_leda_d3_rat_scaled_vector<Self>              Construct_scaled_vector_3;
   typedef Construct_leda_d3_rat_translated_point<Self>           Construct_translated_point_3;
   typedef Construct_leda_d3_rat_point_on<Self>                   Construct_point_on_3;
   typedef Construct_leda_d3_rat_projected_point<Self>            Construct_projected_point_3;
   typedef Construct_leda_d3_rat_lifted_point<Self>               Construct_lifted_point_3;
   typedef Construct_leda_d3_rat_vertex<Self>                     Construct_vertex_3;
   typedef Construct_leda_d3_rat_supporting_line<Self>            Construct_supporting_line_3;
   typedef Construct_leda_d3_rat_supporting_plane<Self>           Construct_supporting_plane_3;
   typedef Construct_leda_d3_rat_orthogonal_vector<Self>          Construct_orthogonal_vector_3;
   typedef Construct_leda_d3_rat_base_vector<Self>                Construct_base_vector_3;
   typedef Construct_leda_d3_rat_perpendicular_plane<Self>        Construct_perpendicular_plane_3;   
   typedef Construct_leda_d3_rat_perpendicular_line<Self>         Construct_perpendicular_line_3;
   typedef Construct_leda_d3_rat_midpoint<Self>                   Construct_midpoint_3;
   typedef Construct_leda_d3_rat_center<Self>                     Construct_center_3;
   typedef Construct_leda_d3_rat_centroid<Self>                   Construct_centroid_3;         
   typedef Construct_leda_d3_rat_circumcenter<Self>               Construct_circumcenter_3;
   typedef Construct_leda_d3_rat_cross_product_vector<Self>       Construct_cross_product_vector_3;
   typedef Construct_leda_d3_rat_opposite_direction<Self>         Construct_opposite_direction_3;
   typedef Construct_leda_d3_rat_opposite_segment<Self>           Construct_opposite_segment_3;
   typedef Construct_leda_d3_rat_opposite_ray<Self>               Construct_opposite_ray_3;
   typedef Construct_leda_d3_rat_opposite_line<Self>              Construct_opposite_line_3;       
   typedef Construct_leda_d3_rat_opposite_plane<Self>             Construct_opposite_plane_3;
   typedef Construct_leda_d3_rat_opposite_sphere<Self>            Construct_opposite_sphere_3;
   typedef Construct_leda_d3_rat_opposite_vector<Self>            Construct_opposite_vector_3;
   
   // moved from 2d kernel ...
   typedef Construct_leda_rat_projected_xy_point<Self>            Construct_projected_xy_point_2;
   
   // 3d assign/intersection
   typedef Assign_leda_rat_3                                Assign_3;
   typedef CGAL_intersect_leda_rat_3<HELP_KERNEL>           Intersect_3;
   
   // 3d computations
   typedef CGAL_compute_leda_d3_rat_squared_distance<HELP_KERNEL>  Compute_squared_distance_3;   
   typedef Compute_leda_d3_rat_squared_length<Self>                Compute_squared_length_3;
   typedef Compute_leda_d3_rat_squared_radius<Self>                Compute_squared_radius_3;
   typedef Compute_leda_d3_rat_squared_area<Self>                  Compute_squared_area_3;  
   typedef Compute_leda_d3_rat_volume<Self>                        Compute_volume_3;

   // 3d generalized predicates ...
   typedef Predicate_leda_d3_rat_angle<Self>                      Angle_3;
   typedef Predicate_leda_d3_rat_equal<Self>                      Equal_3;
   typedef Predicate_leda_d3_rat_equal_x<Self>                    Equal_x_3;
   typedef Predicate_leda_d3_rat_equal_y<Self>                    Equal_y_3;
   typedef Predicate_leda_d3_rat_equal_z<Self>                    Equal_z_3;
   typedef Predicate_leda_d3_rat_equal_xy<Self>                   Equal_xy_3;
   
   // ----------------------------------------------------------------------------------------
   // Equal_xyz_3 is undocumented, but in the interface macros ...
   // ----------------------------------------------------------------------------------------  
   typedef Predicate_leda_d3_rat_equal_xyz<Self>                  Equal_xyz_3;      
   typedef Predicate_leda_d3_rat_less_x<Self>                     Less_x_3;
   typedef Predicate_leda_d3_rat_less_y<Self>                     Less_y_3;
   typedef Predicate_leda_d3_rat_less_z<Self>                     Less_z_3;
   typedef Predicate_leda_d3_rat_less_xy<Self>                    Less_xy_3;
   typedef Predicate_leda_d3_rat_less_xyz<Self>                   Less_xyz_3;   
   typedef Predicate_leda_d3_rat_compare_x<Self>                  Compare_x_3;
   typedef Predicate_leda_d3_rat_compare_y<Self>                  Compare_y_3;
   typedef Predicate_leda_d3_rat_compare_z<Self>                  Compare_z_3; 
   typedef Predicate_leda_d3_rat_compare_xy<Self>                 Compare_xy_3;
   typedef Predicate_leda_d3_rat_compare_xyz<Self>                Compare_xyz_3;   
   typedef Predicate_leda_d3_rat_less_signed_distance_to_plane<Self>  Less_signed_distance_to_plane_3; 
   typedef Predicate_leda_d3_rat_less_distance_to_point<Self>     Less_distance_to_point_3;  
   typedef Predicate_leda_d3_rat_compare_distance<Self>           Compare_distance_3;
   typedef Predicate_leda_d3_rat_collinear<Self>                  Collinear_3;
   typedef Predicate_leda_d3_rat_coplanar<Self>                   Coplanar_3;   
   typedef Predicate_leda_d3_rat_orientation<Self>                Orientation_3;
   typedef Predicate_leda_d3_rat_coplanar_orientation<Self>       Coplanar_orientation_3;   
   typedef Predicate_leda_d3_rat_coplanar_side_of_bounded_circle<Self>   Coplanar_side_of_bounded_circle_3;
   typedef Predicate_leda_d3_rat_side_of_oriented_sphere<Self>    Side_of_oriented_sphere_3; 
   typedef Predicate_leda_d3_rat_side_of_bounded_sphere<Self>     Side_of_bounded_sphere_3;
   typedef Predicate_leda_d3_rat_is_degenerate<Self>              Is_degenerate_3;
   typedef Predicate_leda_d3_rat_has_on<Self>                     Has_on_3;
   typedef Predicate_leda_d3_rat_has_on_bounded_side<Self>        Has_on_bounded_side_3;
   typedef Predicate_leda_d3_rat_has_on_unbounded_side<Self>      Has_on_unbounded_side_3;
   typedef Predicate_leda_d3_rat_has_on_boundary<Self>            Has_on_boundary_3;
   typedef Predicate_leda_d3_rat_has_on_positive_side<Self>       Has_on_positive_side_3;
   typedef Predicate_leda_d3_rat_has_on_negative_side<Self>       Has_on_negative_side_3;
   typedef Predicate_leda_d3_rat_oriented_side<Self>              Oriented_side_3;
   typedef Predicate_leda_d3_rat_bounded_side<Self>               Bounded_side_3;
   typedef Predicate_leda_d3_rat_are_ordered_along_line<Self>     Are_ordered_along_line_3;
   typedef Predicate_leda_d3_rat_are_strictly_ordered_along_line<Self>  Are_strictly_ordered_along_line_3;
   typedef Predicate_leda_d3_rat_collinear_are_ordered_along_line<Self> Collinear_are_ordered_along_line_3;
   typedef Predicate_leda_d3_rat_collinear_are_strictly_ordered_along_line<Self> Collinear_are_strictly_ordered_along_line_3;
   typedef Predicate_leda_d3_rat_do_intersect<Self>               Do_intersect_3;
      
#endif

  // ----------------------------------------------------------------------------------------------------------------
  //  constructors
  // ----------------------------------------------------------------------------------------------------------------

  leda_rat_kernel_traits() { }
  leda_rat_kernel_traits(const leda_rat_kernel_traits& kt) { }


  // ----------------------------------------------------------------------------------------------------------------
  //    access functions
  // ----------------------------------------------------------------------------------------------------------------

  // 2d construction objects ...
    
  Construct_point_2
  construct_point_2_object() const
  { return Construct_point_2(); }  

  Construct_vector_2
  construct_vector_2_object() const
  { return Construct_vector_2(); }   
  
  Construct_direction_2
  construct_direction_2_object() const
  { return Construct_direction_2(); }  

  Construct_direction_of_line_2
  construct_direction_of_line_2_object() const
  { return Construct_direction_of_line_2(); } 
  
  Construct_segment_2
  construct_segment_2_object() const
  { return Construct_segment_2(); }

  Construct_line_2
  construct_line_2_object() const
  { return Construct_line_2(); }
  
  Construct_ray_2
  construct_ray_2_object() const
  { return Construct_ray_2(); }  

  Construct_circle_2
  construct_circle_2_object() const
  { return Construct_circle_2(); } 

  Construct_triangle_2
  construct_triangle_2_object() const
  { return Construct_triangle_2(); } 

  Construct_iso_rectangle_2
  construct_iso_rectangle_2_object() const
  { return Construct_iso_rectangle_2(); } 
  
  Construct_object_2
  construct_object_2_object() const
  { return Construct_object_2(); }   
  
  Construct_scaled_vector_2
  construct_scaled_vector_2_object() const
  { return Construct_scaled_vector_2(); }   

  Construct_translated_point_2
  construct_translated_point_2_object() const
  { return Construct_translated_point_2(); } 

  Construct_point_on_2
  construct_point_on_2_object() const
  { return Construct_point_on_2(); } 
  
  Construct_projected_point_2
  construct_projected_point_2_object() const
  { return Construct_projected_point_2(); }   
  
  // Construct_projected_xy_point_2 goes
  // to the 3d constructions ...
    
  Construct_vertex_2
  construct_vertex_2_object() const
  { return  Construct_vertex_2(); }

  Construct_supporting_line_2
  construct_supporting_line_2_object() const
  { return Construct_supporting_line_2(); } 
  
  Construct_perpendicular_vector_2
  construct_perpendicular_vector_2_object() const
  { return Construct_perpendicular_vector_2(); }   
  
  Construct_perpendicular_direction_2
  construct_perpendicular_direction_2_object() const
  { return Construct_perpendicular_direction_2(); }   
  
  Construct_perpendicular_line_2
  construct_perpendicular_line_2_object() const
  { return Construct_perpendicular_line_2(); }   

  Construct_midpoint_2
  construct_midpoint_2_object() const
  { return Construct_midpoint_2(); }   

  Construct_center_2
  construct_center_2_object() const
  { return Construct_center_2(); }
  
  Construct_centroid_2
  construct_centroid_2_object() const
  { return Construct_centroid_2(); }  
    
  Construct_circumcenter_2
  construct_circumcenter_2_object() const
  { return Construct_circumcenter_2(); } 
  
  Construct_bisector_2
  construct_bisector_2_object() const
  { return Construct_bisector_2(); }    

  Construct_opposite_direction_2
  construct_opposite_direction_2_object() const
  { return Construct_opposite_direction_2(); }  
  
  Construct_opposite_segment_2
  construct_opposite_segment_2_object() const
  { return Construct_opposite_segment_2(); }    

  Construct_opposite_ray_2
  construct_opposite_ray_2_object() const
  { return Construct_opposite_ray_2(); }  

  Construct_opposite_line_2
  construct_opposite_line_2_object() const
  { return Construct_opposite_line_2(); }  

  Construct_opposite_triangle_2
  construct_opposite_triangle_2_object() const
  { return Construct_opposite_triangle_2(); }  
  
  Construct_opposite_circle_2
  construct_opposite_circle_2_object() const
  { return Construct_opposite_circle_2(); }  
  
  Construct_opposite_vector_2
  construct_opposite_vector_2_object() const
  { return Construct_opposite_vector_2(); }      

  // ----------------------------------------------------------------------
  // 2d intersections and assignments
  // ----------------------------------------------------------------------

  Assign_2 
  assign_2_object() const
  { return  Assign_2(); }
  
  Intersect_2
  intersect_2_object() const
  { return   Intersect_2(); }
  

  // ----------------------------------------------------------------------
  // 2d computations
  // ----------------------------------------------------------------------

  Compute_squared_distance_2
  compute_squared_distance_2_object() const
  { return Compute_squared_distance_2(); }
  
  Compute_squared_length_2
  compute_squared_length_2_object() const
  { return Compute_squared_length_2(); }
  
  Compute_squared_radius_2
  compute_squared_radius_2_object() const
  { return Compute_squared_radius_2(); }
  
  Compute_area_2 
  compute_area_2_object() const
  { return Compute_area_2(); }

  Compute_y_at_x_2
  compute_y_at_x_2_object() const
  { return Compute_y_at_x_2(); }

  // ----------------------------------------------------------------------
  // 2d predicates
  // ----------------------------------------------------------------------

  Angle_2
  angle_2_object() const
  { return Angle_2(); }

  Equal_2
  equal_2_object() const 
  { return Equal_2(); } 

  Equal_x_2
  equal_x_2_object() const 
  { return Equal_x_2(); } 
  
  Equal_y_2
  equal_y_2_object() const 
  { return Equal_y_2(); }   
  
  Less_x_2
  less_x_2_object() const 
  { return Less_x_2(); }   

  Less_y_2
  less_y_2_object() const 
  { return Less_y_2(); }  
  
  Less_xy_2
  less_xy_2_object() const 
  { return Less_xy_2(); } 

  Less_yx_2
  less_yx_2_object() const 
  { return Less_yx_2(); } 

  Compare_x_2
  compare_x_2_object() const 
  { return Compare_x_2(); } 

  Compare_x_at_y_2
  compare_x_at_y_2_object() const 
  { return Compare_x_at_y_2(); } 
  
  Compare_y_2
  compare_y_2_object() const 
  { return Compare_y_2(); } 

  Compare_xy_2
  compare_xy_2_object() const 
  { return Compare_xy_2(); } 

  Compare_y_at_x_2
  compare_y_at_x_2_object() const 
  { return Compare_y_at_x_2(); } 

  Compare_distance_2
  compare_distance_2_object() const 
  { return Compare_distance_2(); } 

  Compare_angle_with_x_axis_2
  compare_angle_with_x_axis_2_object() const 
  { return Compare_angle_with_x_axis_2(); } 
  
  Compare_slope_2
  compare_slope_2_object() const
  { return Compare_slope_2(); }

  Less_distance_to_point_2
  less_distance_to_point_2_object( ) const
  { return Less_distance_to_point_2( ); } 

  Less_signed_distance_to_line_2
  less_signed_distance_to_line_2_object( ) const
  { return Less_signed_distance_to_line_2( ); } 

  Less_rotate_ccw_2
  less_rotate_ccw_2_object( ) const
  { return Less_rotate_ccw_2( ); }

  Leftturn_2
  leftturn_2_object() const
  { return Leftturn_2(); }
  
  Left_turn_2
  left_turn_2_object() const
  { return Left_turn_2(); }  

  Collinear_2
  collinear_2_object() const
  { return Collinear_2(); }
  
  Orientation_2
  orientation_2_object() const
  { return Orientation_2(); } 
  
  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  { return Side_of_oriented_circle_2(); } 
  
  Side_of_bounded_circle_2
  side_of_bounded_circle_2_object() const
  { return Side_of_bounded_circle_2(); }
  
  Is_horizontal_2
  is_horizontal_2_object() const
  { return Is_horizontal_2(); }
  
  Is_vertical_2
  is_vertical_2_object() const
  { return Is_vertical_2(); }  
    
  Is_degenerate_2
  is_degenerate_2_object() const
  { return Is_degenerate_2(); } 
  
  Has_on_2
  has_on_2_object() const
  { return Has_on_2(); }    
  
  Collinear_has_on_2
  collinear_has_on_2_object() const
  { return Collinear_has_on_2(); } 
  
  Has_on_bounded_side_2
  has_on_bounded_side_2_object() const
  { return Has_on_bounded_side_2(); }
  
  Has_on_unbounded_side_2
  has_on_unbounded_side_2_object() const
  { return Has_on_unbounded_side_2(); }       

  Has_on_boundary_2
  has_on_boundary_2_object() const
  { return Has_on_boundary_2(); }
  
  Has_on_positive_side_2
  has_on_positive_side_2_object() const
  { return Has_on_positive_side_2(); }  

  Has_on_negative_side_2
  has_on_negative_side_2_object() const
  { return Has_on_negative_side_2(); }  

  Oriented_side_2
  oriented_side_2_object() const
  { return Oriented_side_2(); }  
  
  Bounded_side_2
  bounded_side_2_object() const
  { return Bounded_side_2(); }    

  Are_ordered_along_line_2
  are_ordered_along_line_2_object() const
  { return Are_ordered_along_line_2(); }  
  
  Are_strictly_ordered_along_line_2
  are_strictly_ordered_along_line_2_object() const
  { return Are_strictly_ordered_along_line_2(); }    

  Collinear_are_ordered_along_line_2
  collinear_are_ordered_along_line_2_object() const
  { return Collinear_are_ordered_along_line_2(); }  

  Collinear_are_strictly_ordered_along_line_2
  collinear_are_strictly_ordered_along_line_2_object() const
  { return Collinear_are_strictly_ordered_along_line_2(); }  

  Counterclockwise_in_between_2
  counterclockwise_in_between_2_object() const
  { return Counterclockwise_in_between_2(); }  

  Do_intersect_2
  do_intersect_2_object() const
  { return Do_intersect_2(); } 
  
  // extensions ...
  Do_intersect_to_right_2
  do_intersect_to_right_2_object() const
  { return Do_intersect_to_right_2(); }   
  
  Do_intersect_to_left_2
  do_intersect_to_left_2_object() const
  { return Do_intersect_to_left_2(); }   
  
  
// -----------------------------------------------------------
// support for 3d kernel traits ...
#if defined(CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3)

   // 3d constructions

   Construct_point_3
   construct_point_3_object() const
   { return Construct_point_3(); }
   
   Construct_vector_3
   construct_vector_3_object() const
   { return Construct_vector_3(); }   
   
   Construct_direction_3
   construct_direction_3_object() const
   { return Construct_direction_3(); }   
   
   Construct_direction_of_line_3
   construct_direction_of_line_3_object() const
   { return Construct_direction_of_line_3(); }
   
   Construct_plane_3
   construct_plane_3_object() const
   { return Construct_plane_3(); }
   
   Construct_iso_cuboid_3
   construct_iso_cuboid_3_object() const
   { return Construct_iso_cuboid_3(); }
   
   Construct_line_3
   construct_line_3_object() const
   { return Construct_line_3(); }   
   
   Construct_ray_3
   construct_ray_3_object() const
   { return Construct_ray_3(); } 
   
   Construct_sphere_3
   construct_sphere_3_object() const
   { return Construct_sphere_3(); } 
         
   Construct_segment_3
   construct_segment_3_object() const
   { return Construct_segment_3(); }   
   
   Construct_tetrahedron_3
   construct_tetrahedron_3_object() const
   { return Construct_tetrahedron_3(); }   

   Construct_triangle_3
   construct_triangle_3_object() const
   { return Construct_triangle_3(); }  
   
   Construct_object_3
   construct_object_3_object() const
   { return Construct_object_3(); }   
   
   Construct_scaled_vector_3
   construct_scaled_vector_3_object() const
   { return Construct_scaled_vector_3(); } 
   
   Construct_translated_point_3
   construct_translated_point_3_object() const
   { return Construct_translated_point_3(); }   
   
   Construct_point_on_3
   construct_point_on_3_object() const
   { return Construct_point_on_3(); }    
   
   Construct_projected_point_3
   construct_projected_point_3_object() const
   { return Construct_projected_point_3(); }       
   
   Construct_lifted_point_3
   construct_lifted_point_3_object() const
   { return Construct_lifted_point_3(); }            

   Construct_vertex_3
   construct_vertex_3_object() const
   { return Construct_vertex_3(); } 	 
	 
   Construct_supporting_line_3
   construct_supporting_line_3_object() const
   { return Construct_supporting_line_3(); } 

   Construct_supporting_plane_3
   construct_supporting_plane_3_object() const
   { return Construct_supporting_plane_3(); } 
   
   Construct_orthogonal_vector_3
   construct_orthogonal_vector_3_object() const
   { return Construct_orthogonal_vector_3(); }    

   Construct_base_vector_3
   construct_base_vector_3_object() const
   { return Construct_base_vector_3(); }    

   Construct_perpendicular_plane_3
   construct_perpendicular_plane_3_object() const
   { return Construct_perpendicular_plane_3(); } 
   
   Construct_perpendicular_line_3
   construct_perpendicular_line_3_object() const
   { return Construct_perpendicular_line_3(); }   
 
   Construct_midpoint_3
   construct_midpoint_3_object() const
   { return Construct_midpoint_3(); }  
 
   Construct_center_3
   construct_center_3_object() const
   { return Construct_center_3(); }
   
   Construct_centroid_3
   construct_centroid_3_object() const
   { return Construct_centroid_3(); }   
      
   Construct_circumcenter_3
   construct_circumcenter_3_object() const
   { return Construct_circumcenter_3(); }   
   
   Construct_cross_product_vector_3
   construct_cross_product_vector_3_object() const
   { return Construct_cross_product_vector_3(); }   
      
   Construct_opposite_direction_3
   construct_opposite_direction_3_object() const
   { return Construct_opposite_direction_3(); }   
   
   Construct_opposite_segment_3
   construct_opposite_segment_3_object() const
   { return Construct_opposite_segment_3(); }   

   Construct_opposite_ray_3
   construct_opposite_ray_3_object() const
   { return Construct_opposite_ray_3(); }   

   Construct_opposite_line_3
   construct_opposite_line_3_object() const
   { return Construct_opposite_line_3(); }   

   Construct_opposite_plane_3
   construct_opposite_plane_3_object() const
   { return Construct_opposite_plane_3(); }
   
   Construct_opposite_sphere_3
   construct_opposite_sphere_3_object() const
   { return Construct_opposite_sphere_3(); }   

   Construct_opposite_vector_3
   construct_opposite_vector_3_object() const
   { return Construct_opposite_vector_3(); }   

   // moved from 2d ...
   Construct_projected_xy_point_2
   construct_projected_xy_point_2_object() const           
   { return Construct_projected_xy_point_2(); }      
	          
   // 3d intersections
   
   Assign_3 
   assign_3_object() const
   { return  Assign_3(); }
  
   Intersect_3
   intersect_3_object() const
   { return   Intersect_3(); }   
   
   
   // 3d computations

   Compute_squared_distance_3
   compute_squared_distance_3_object() const
   { return Compute_squared_distance_3(); }
  
   Compute_squared_length_3
   compute_squared_length_3_object() const
   { return Compute_squared_length_3(); }
  
   Compute_squared_radius_3
   compute_squared_radius_3_object() const
   { return Compute_squared_radius_3(); }
  
   Compute_squared_area_3 
   compute_squared_area_3_object() const
   { return Compute_squared_area_3(); }   

   Compute_volume_3 
   compute_volume_3_object() const
   { return Compute_volume_3(); }   
   
   // 3d generalized predicates 
   
   Angle_3
   angle_3_object() const 
   { return Angle_3(); }   
   
   Equal_3
   equal_3_object() const 
   { return Equal_3(); }
   
   Equal_x_3
   equal_x_3_object() const 
   { return Equal_x_3(); }   
   
   Equal_y_3
   equal_y_3_object() const 
   { return Equal_y_3(); }
   
   Equal_z_3
   equal_z_3_object() const 
   { return Equal_z_3(); }
   
   Equal_xy_3
   equal_xy_3_object() const 
   { return Equal_xy_3(); }  
   
   Less_x_3
   less_x_3_object() const 
   { return Less_x_3(); }                  

   Less_y_3
   less_y_3_object() const 
   { return Less_y_3(); }  
   
   Less_z_3
   less_z_3_object() const 
   { return Less_z_3(); }     
   
   Less_xy_3
   less_xy_3_object() const 
   { return Less_xy_3(); }     
   
   Less_xyz_3
   less_xyz_3_object() const 
   { return Less_xyz_3(); }     
     
   Compare_x_3
   compare_x_3_object() const 
   { return Compare_x_3(); }   
   
   Compare_y_3
   compare_y_3_object() const 
   { return Compare_y_3(); }      
   
   Compare_z_3 
   compare_z_3_object() const 
   { return Compare_z_3(); }  
   
   Compare_xy_3
   compare_xy_3_object() const 
   { return Compare_xy_3(); }   
   
   Compare_xyz_3
   compare_xyz_3_object() const 
   { return Compare_xyz_3(); }          

   Less_signed_distance_to_plane_3
   less_signed_distance_to_plane_3_object() const
   { return Less_signed_distance_to_plane_3(); }
   
   Less_distance_to_point_3
   less_distance_to_point_3_object() const
   { return Less_distance_to_point_3(); }
   
   Compare_distance_3
   compare_distance_3_object() const 
   { return Compare_distance_3(); }      

   Collinear_3
   collinear_3_object() const 
   { return Collinear_3(); }    
   
   Coplanar_3
   coplanar_3_object() const 
   { return Coplanar_3(); }       
   
   Orientation_3
   orientation_3_object() const 
   { return Orientation_3(); }      
   
   Coplanar_orientation_3
   coplanar_orientation_3_object() const 
   { return Coplanar_orientation_3(); }   
   
   Coplanar_side_of_bounded_circle_3
   coplanar_side_of_bounded_circle_3_object() const 
   { return Coplanar_side_of_bounded_circle_3(); }
      
   Side_of_oriented_sphere_3
   side_of_oriented_sphere_3_object() const 
   { return Side_of_oriented_sphere_3(); }   
 
   Side_of_bounded_sphere_3
   side_of_bounded_sphere_3_object() const 
   { return Side_of_bounded_sphere_3(); }    
   
   Is_degenerate_3
   is_degenerate_3_object() const
   { return Is_degenerate_3(); }
   
   Has_on_3
   has_on_3_object() const
   { return Has_on_3(); }
   
   Has_on_bounded_side_3
   has_on_bounded_side_3_object() const
   { return Has_on_bounded_side_3(); }
   
   Has_on_unbounded_side_3
   has_on_unbounded_side_3_object() const
   { return Has_on_unbounded_side_3(); }
   
   Has_on_boundary_3
   has_on_boundary_3_object() const
   { return Has_on_boundary_3(); }
   
   Has_on_positive_side_3
   has_on_positive_side_3_object() const
   { return Has_on_positive_side_3(); }
   
   Has_on_negative_side_3
   has_on_negative_side_3_object() const
   { return Has_on_negative_side_3(); }
   
   Oriented_side_3
   oriented_side_3_object() const
   { return Oriented_side_3(); }
   
   Bounded_side_3
   bounded_side_3_object() const
   { return Bounded_side_3(); }
   
   Are_ordered_along_line_3
   are_ordered_along_line_3_object() const
   { return Are_ordered_along_line_3(); }
   
   Are_strictly_ordered_along_line_3
   are_strictly_ordered_along_line_3_object() const
   { return Are_strictly_ordered_along_line_3(); }
   
   Collinear_are_ordered_along_line_3
   collinear_are_ordered_along_line_3_object() const
   { return Collinear_are_ordered_along_line_3(); }
   
   Collinear_are_strictly_ordered_along_line_3
   collinear_are_strictly_ordered_along_line_3_object() const
   { return Collinear_are_strictly_ordered_along_line_3(); }
   
   Do_intersect_3 
   do_intersect_3_object() const  
   { return Do_intersect_3(); } 
   
#endif

};

CGAL_END_NAMESPACE

#endif


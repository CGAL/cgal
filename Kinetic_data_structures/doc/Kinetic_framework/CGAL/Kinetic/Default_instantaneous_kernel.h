
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsFrameworkClasses

This class provides a model of the `Kinetic::InstantaneousKernel` for 
use with general %Cartesian Geometry. It provides all the predicates 
needed for Delaunay triangulations and regular triangulations. 

\cgalModels `Kinetic::InstantaneousKernel`

*/
template< typename SimulationTraits >
class Default_instantaneous_kernel {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef unspecified_type Orientation_2; 

/*!

*/ 
typedef unspecified_type Orientation_3; 

/*!

*/ 
typedef unspecified_type Side_of_oriented_circle_2; 

/*!

*/ 
typedef unspecified_type Side_of_oriented_sphere_3; 

/*!

*/ 
typedef unspecified_type Power_test_3; 

/*!

*/ 
typedef unspecified_type Weighted_orientation_3; 

/*!

*/ 
typedef unspecified_type Compare_x_1; 

/*!

*/ 
typedef unspecified_type Compare_x_2; 

/*!

*/ 
typedef unspecified_type Compare_y_2; 

/*!

*/ 
typedef unspecified_type Compare_x_3; 

/*!

*/ 
typedef unspecified_type Compare_y_3; 

/*!

*/ 
typedef unspecified_type Compare_z_3; 

/*!

*/ 
typedef unspecified_type Compare_distance_2; 

/*!

*/ 
typedef unspecified_type Compare_distance_3; 

/*!
Note that this one does not work if the current time is not a `NT`. 
*/ 
typedef unspecified_type Coplanar_orientation_3; 

/*!
Note that this one does not work if the current time is not a `NT`. 
*/ 
typedef unspecified_type Coplanar_side_of_bounded_circle_3; 

/// @}

}; /* end Kinetic::Default_instantaneous_kernel */
} /* end namespace Kinetic */
} /* end namespace CGAL */

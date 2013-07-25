
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsFrameworkClasses

This class provides a model of `Kinetic::Kernel` for use with general %Cartesian geometry. 

The IO format for points is currently \f$ p_0\f$, \f$ p_1\f$, ... \f$ w\f$. \f$ p_i\f$ and \f$ w\f$ are instances of Function. There IO format is typically \f$ c_0+c_1t+c_2t^2+...\f$. Beware of issues with \cgal IO of the coeffients as exact number typles often require that the coefficents be expressed as \f$ a/b\f$ even when \f$ b\f$ is 1. 

\cgalModels `Kinetic::Kernel`

*/
template< typename FunctionKernel >
class Cartesian {
public:

/// \name Types 
/// @{

/*!
This is a model of `Kinetic::Certificate`. 
*/ 
typedef unspecified_type Certificate; 

/*!

*/ 
typedef unspecified_type Point_1; 

/*!

*/ 
typedef unspecified_type Point_2; 

/*!

*/ 
typedef unspecified_type Point_3; 

/*!

*/ 
typedef unspecified_type Weighted_point_3; 

/// @}

/// \name Certificate Functors
/// The following are functors which generate `Certificate`
/// objects. Each has a corresponding `_object` method which creates
/// the functor. They are models of `Kinetic::CertificateGenerator`.
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

/// @}

}; /* end Kinetic::Cartesian */
} /* end namespace Kinetic */
} /* end namespace CGAL */

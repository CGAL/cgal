
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsSupport

The class `Kinetic::Enclosing_box_2` keeps the points in the simulation inside of a 
box. Whenever the points come close to the wall of the box they bounce off of the wall. 

Note that, in general, points hit the wall of the box at times which 
are not easily represented by standard (rational) number types. The 
resulting trajectories would also have non-rational coefficients, 
complicating and slowing the simulation. In order to handle this, the 
`Kinetic::Enclosing_box_2` bounces the points at the nearest easily representable 
time before the point would leave the box. 

*/
template< typename Traits >
class Enclosing_box_2 {
public:

/// \name Types 
/// @{

/*!
The number type used to represent the walls of the box and perform calculations. Generally this is `Traits::NT`. 
*/ 
typedef unspecified_type NT; 

/// @} 

/// \name Creation 
/// @{

/*!
This constructs a bounding box with the dimensions specified by the last 4 arguments. They are optional and will take the values \f$ \pm\f$ 10 if omitted. 
*/ 
Enclosing_box_2(Traits, NT xmin, NT xmax, NT ymin, NT ymax); 

/// @}

}; /* end Kinetic::Enclosing_box_2 */
} /* end namespace Kinetic */
} /* end namespace CGAL */

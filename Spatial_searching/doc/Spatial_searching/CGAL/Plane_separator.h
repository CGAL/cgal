namespace CGAL {

/*!
\ingroup AdvancedClasses

The class `Plane_separator` implements a plane separator, i.e., a
hyperplane that is used to separate two half spaces. \advanced This
hyperplane is defined by a cutting dimension \f$ d\f$ and a cutting
value \f$ v\f$ as \f$ x_d=v\f$, where \f$ v\f$ denotes the \f$
d^{th}\f$ coordinate value.

\models ::SpatialSeparator 

*/
template< typename FT >
class Plane_separator {
public:

/// \name Creation 
/// @{

/*! 
Constructs a separator that separates two half spaces by a hyperplane 
defined by \f$ x_d=v\f$, where \f$ v\f$ denotes the \f$ d^{th}\f$ coordinate value. 
*/ 
Plane_separator(int d, FT v); 

/*! 
Copy constructor. 
*/ 
Plane_separator(Plane_separator<FT> p); 

/// @} 

/// \name Operations 
/// @{

/*! 
Sets the cutting dimension to `d`. 
*/ 
void set_cutting_dimension(int d); 

/*! 
Sets the cutting value to `v`. 
*/ 
void set_cutting_value(FT v); 

/*! 
Returns the number of the cutting dimension. 
*/ 
int cutting_dimension() const; 

/*! 
Returns the cutting value. 
*/ 
FT cutting_value() const; 

/*! 
Assignment operator. 
*/ 
Plane_separator<FT> operator=(Plane_separator<FT> s2); 

/// @}

}; /* end Plane_separator */

/*! 
Inserts the plane separator `s` in the output stream `os` and returns `os`. 
\relates Plane_separator 
*/ 
template<class FT> 
std::ostream& operator<<(std::ostream& os, Plane_separator<FT> s); 

} /* end namespace CGAL */

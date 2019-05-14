namespace CGAL {

/*!
\ingroup PkgPlacementOfStreamlines2


This class provides a 2D vector field specified by a set of sample 
points defined on a regular grid, with a bilinear interpolation scheme 
over its cells (i.e.\ for each point `p` in a cell `c`, the 
vector value is interpolated from the vertices of `c`). 

\tparam StreamLinesTraits_2 has to instantiated by a model of the concept `StreamLinesTraits_2`. 

\cgalModels `VectorField_2`

\sa `Triangular_field_2<StreamLinesTraits_2>` 

*/
template< typename StreamLinesTraits_2 >
class Regular_grid_2 {
public:

/// \name Types 
/// @{

/*!
the scalar type. 
*/ 
typedef StreamLinesTraits_2::FT FT; 

/*!
the point type. 
*/ 
typedef StreamLinesTraits_2::Point_2 Point_2; 

/*!
the vector type. 
*/ 
typedef StreamLinesTraits_2::Vector_2 Vector_2; 

/// @} 

/// \name Creation 
/// @{

/*!
Generate a regular grid  whose size is `x_size` by `y_size`, while `x_samples` and 
`y_samples` specify the number of samples on `x` and `y`. 
*/ 
Regular_grid_2(int x_samples, int y_samples, FT x_size, FT y_size); 

/// @} 

/// \name Modifiers 
/// In addition to the minimum interface required by the concept
/// definition, the class `Regular_grid_2` provides the following
/// function to fill the vector field with the user data.
/// @{

/*!
Attribute the vector v 
to the position (i,j) on the regular grid. 
*/ 
void set_xy(int i, int j, Vector_2 v); 

/// @} 

/// \name Access Functions 
/// @{

/*!
returns the dimension of the grid. 
*/ 
std::pair<int, int> get_dimension(); 

/*!
returns the size of the grid. 
*/ 
std::pair<FT, FT> get_size(); 

/// @}

}; /* end Regular_grid_2 */
} /* end namespace CGAL */

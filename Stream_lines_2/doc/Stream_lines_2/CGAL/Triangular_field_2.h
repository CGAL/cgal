namespace CGAL {

/*!
\ingroup PkgPlacementOfStreamlines2


This class provides a vector field specified by a set of sample points 
defined on a triangulated domain. All sample points are inserted to a 
`Delaunay triangulation`, and for each point `p` in the domain 
located in a face `f`, its vector value is interpolated from the 
vertices of the face `f`. 

\tparam StreamLinesTraits_2 has to be instantiated by a model of the concept `StreamLinesTraits_2`. 

\cgalModels `VectorField_2`

\sa `Regular_grid_2<StreamLinesTraits_2>` 

*/
template< typename StreamLinesTraits_2 >
class Triangular_field_2 {
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
Defines the points in the range 
`[first_point, `last_point) 
as the sample points of the grid, with the corresponding number of vectors started at `first_vector`. 
\pre The `value_type` of `InputIterator1` is `Point_2`. 
\pre The `value_type` of `InputIterator2` is `Vector_2`. 
*/ 
Triangular_field_2( InputIterator1 first_point, InputIterator1 last_point, InputIterator2 first_vector); 

/// @}

}; /* end Triangular_field_2 */
} /* end namespace CGAL */

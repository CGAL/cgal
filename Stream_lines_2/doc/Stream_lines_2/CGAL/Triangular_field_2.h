namespace CGAL {

/*!
\ingroup PkgStreamLines2Ref


This class provides a vector field specified by a set of sample points
defined on a triangulated domain. All sample points are inserted to a
Delaunay triangulation, and for each point `p` in the domain
located in a face `f`, its vector value is interpolated from the
vertices of the face `f`.

\tparam StreamLinesTraits_2 has to be instantiated by a model of the concept `StreamLinesTraits_2`.

\cgalModels{VectorField_2}

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
Defines the points in the range `[first_point, last_point)`
as the sample points of the triangular field, with the corresponding number of vectors started at `first_vector`.
\tparam PointInputIterator must be an input iterator with the value type `Point_2`.
\tparam VectorInputIterator must be an input iterator with the value type `Vector_2`.
*/
mplate <class PointIterator1, class VectorInputIterator>
Triangular_field_2( PointInputIterator first_point, PointInputIterator last_point, VectorInputIterator first_vector);

/// @}

}; /* end Triangular_field_2 */
} /* end namespace CGAL */

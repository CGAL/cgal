namespace CGAL {

/*!
\ingroup PkgStreamLines2Ref

The class `Stream_lines_2` generates a placement of streamlines
in a 2D domain according to a bidimensional vector field.

The class places streamlines according to
a specified density and gives access to the generated streamlines via two
iterators over a container of iterators that provide access to the
streamline points.


\tparam VectorField_2 must be a model of the concept
`VectorField_2`.
\tparam  Integrator_2 is a function object and must be a model of the concept `Integrator_2`.

*/
template< typename VectorField_2, typename Integrator_2 >
class Stream_lines_2 {
public:

/// \name Types
/// @{

/*!
the traits class.
*/
typedef VectorField_2::Geom_traits Geom_traits;

/*!
the scalar type.
*/
typedef VectorField_2::FT FT;

/*!
the point type.
*/
typedef VectorField_2::Point_2 Point_2;

/*!
the vector type.
*/
typedef VectorField_2::Vector_2 Vector_2;

/// @}

/// \name Streamline Iterators
/// The following iterators allow to visit all the streamlines generated
/// by the constructor or the update function.

/// @{

/*!
iterator of points with value type `Point_2`.
*/
typedef unspecified_type Point_iterator_2;

/*!
an iterator to visit the streamlines with value type `std::pair<Point_iterator_2, Point_iterator_2>`.
*/
typedef unspecified_type Stream_line_iterator_2;

/// @}


/*!
Constructor which generates a streamline placement.
*/
Stream_lines_2(VectorField_2 vector_field_2, Integrator_2 integrator_2, FT
separating_distance, FT saturation_ratio);


/// \name Modifiers
/// @{

/*!
Modify the
separating distance.
*/
void set_separating_distance(FT new_value);

/*!
Modify the
saturation ratio.
*/
void set_saturation_ratio(FT new_value);

/*!
Update the placement after changing the
separating distance or the saturation ratio.
*/
void update();

/// @}

/// \name Access Functions
/// @{

/*!
returns the separating distance.
*/
FT get_separating_distance() const;

/*!
returns the saturation ratio.
*/
FT get_saturation_ratio() const;

/*!
prints the streamlines to an ASCII file: line by line, and point by point.
*/
void print_stream_lines(std::ofstream & fw);

/// @}

/// \name Streamline Iterators
/// @{

/*!
Starts at the first streamline
*/
Stream_line_iterator begin() const;

/*!
Past-the-end iterator
*/
Stream_line_iterator end() const;

/// @}

}; /* end Stream_lines_2 */
} /* end namespace CGAL */

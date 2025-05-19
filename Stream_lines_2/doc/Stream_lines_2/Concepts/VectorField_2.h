
/*!
\ingroup PkgStreamLines2Concepts
\cgalConcept

The concept `VectorField_2` describes the set of requirements for
the first template parameter of the class
`CGAL::Stream_lines_2<VectorField_2,Integrator_2>`. This concept
provides the types of the geometric primitives used in the placement
of streamlines and some functions for answering different queries.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Regular_grid_2<StreamLinesTraits_2>}
\cgalHasModels{CGAL::Triangular_field_2<StreamLinesTraits_2>}
\cgalHasModelsEnd

*/

class VectorField_2 {
public:

/// \name Types
/// @{

/*!
The traits class.
*/
typedef unspecified_type Geom_traits;

/*!
The scalar type.
*/
typedef unspecified_type FT;

/*!
The point type.
*/
typedef unspecified_type Point_2;

/*!
The vector type.
*/
typedef unspecified_type Vector_2;

/// @}

/// \name Creation
/// @{

/*!
Any constructor has to allow the user to fill the vector values, i.e., assign a vector to each position within the domain.
*/
VectorField_2();

/// @}

/// \name Query Functions
/// @{

/*!
returns the bounding box of the whole domain.
*/
Geom_traits::Iso_rectangle_2 bbox();

/*!
returns the vector field value and the local density.
\pre `is_in_domain(p)` must be `true`.
*/
std::pair<Vector_2,FT> get_field(Point_2 p);

/*!
returns `true` if the point `p`
is inside the domain boundaries, `false` otherwise.
*/
bool is_in_domain(Point_2 p);

/*!
returns the integration step at the point `p`, i.e., the distance between `p` and the next point in the polyline.
\pre `is_in_domain(p)` must be `true`.
*/
FT get_integration_step(Point_2 p);

/// @}

}; /* end VectorField_2 */


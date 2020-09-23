/*!
\ingroup PkgBarycentricCoordinates2Concepts
\cgalConcept

Requirements of the template parameter `Traits` for all the classes with two-dimensional barycentric coordinates from the namespace `CGAL::Barycentric_coordinates`.

\cgalHasModel All models of `Kernel`

*/

class BarycentricTraits_2 {

public:

/// \name Types
/// @{

/*!
        A model of `FieldNumberType`.
*/
typedef unspecified_type FT;

/// @}

/// \name Two-dimensional Geometric Objects
/// @{

/*!
        A model of `Kernel::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
        A model of `Kernel::Vector_2`.
*/
typedef unspecified_type Vector_2;

/// @}

/// \name Two-dimensional Constructions
/// @{

/*!
        A model of this concept must provide an FT operator(const Point_2 &p, const Point_2 &q, const Point_2 &r)
        that returns the signed area of the triangle defined by the points p, q, and r.
*/
typedef unspecified_type Compute_area_2;

/*!
    A model of this concept must provide an FT operator(const Point_2 &p, const Point_2 &q)
    that returns the squared Euclidean distance between the points p and q.
*/
typedef unspecified_type Compute_squared_distance_2;

/*!
    A model of this concept must provide an FT operator(const Vector_2 &p)
    that returns the squared length of the vector p.
*/
typedef unspecified_type Compute_squared_length_2;

/*!
    A model of this concept must provide an FT operator(const Vector_2 &p, const Vector_2 &q)
    that returns the scalar product between the vectors p and q.
*/
typedef unspecified_type Compute_scalar_product_2;

/// @}

/// \name Two-dimensional Generalized Predicates
/// @{

/*!
        A model of this concept must provide a bool operator(const Point_2 &p, const Point_2 &q)
        that returns true if p = q and false otherwise.
*/
typedef unspecified_type Equal_2;

/*!
        A model of this concept must provide a bool operator(const Point_2 &p, const Point_2 &q, const Point_2 &r)
        that returns true if the points p, q, and r are collinear and false otherwise.
*/
typedef unspecified_type Collinear_2;

/*!
        A model of this concept must provide a bool operator(const Point_2 &p, const Point_2 &q, const Point_2 &r)
        that returns true if the point q lies between the points p and r and all three points are collinear.
*/
typedef unspecified_type Collinear_are_ordered_along_line_2;

/// @}

}; /* end BarycentricTraits_2 */


namespace CGAL {

/*!
\ingroup PkgPolytopeDistanceDRef

The class `Polytope_distance_d_traits_2` is a traits class for the \f$ d\f$-dimensional
optimization algorithms using the two-dimensional \cgal kernel.


\tparam K must be a model for `Kernel`.

\tparam ET must be models for `RingNumberType`. The default is  `K::RT`.
\tparam NT must be models for `RingNumberType`. The default is  `K::RT`.

\cgalModels{PolytopeDistanceDTraits}

\sa `CGAL::Polytope_distance_d<Traits>`
\sa `CGAL::Polytope_distance_d_traits_3<K,ET,NT>`
\sa `CGAL::Polytope_distance_d_traits_d<K,ET,NT>`
\sa `PolytopeDistanceDTraits`

*/
template< typename K, typename ET, typename NT >
class Polytope_distance_d_traits_2 {
public:

/// \name Types
/// @{

/*!
the point type.
*/
typedef K::Point_2 Point_d;

/*!
typedef to `K::Rep_tag`.
*/
typedef unspecified_type Rep_tag;

/*!
the ring type.
*/
typedef K::RT RT;

/*!
the field type.
*/
typedef K::FT FT;

/*!
functor returning `2`.
*/
typedef unspecified_type Access_dimension_d;

/*!
functor constructing the begin iterator of the homogeneous coordinates of a point.
*/
typedef unspecified_type Access_coordinates_begin_d;

/*!
functor constructing a point from a coordinate range.
*/
typedef unspecified_type Construct_point_d;

/*!
second template parameter (default is `K::RT`).
*/
typedef unspecified_type ET;

/*!
third template parameter (default is `K::RT`).
*/
typedef unspecified_type NT;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Polytope_distance_d_traits_2( );

/*!
copy constructor.
*/
Polytope_distance_d_traits_2(
const Polytope_distance_d_traits_2<K,ET,NT>&);

/// @}

/// \name Operations
/// The following functions just return the corresponding function class object.
/// @{

/*!

*/
Access_dimension_d
access_dimension_d_object() const;

/*!

*/
Access_coordinates_begin_d
access_coordinates_begin_d_object() const;

/*!

*/
Construct_point_d
construct_point_d_object() const;

/// @}

}; /* end Polytope_distance_d_traits_2 */
} /* end namespace CGAL */

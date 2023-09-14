
namespace CGAL {

/*!
\ingroup PkgPolytopeDistanceDRef

The class `Polytope_distance_d_traits_d` is a traits class for the \f$ d\f$-dimensional
optimisation algorithms using the \f$ d\f$-dimensional \cgal kernel.


\tparam K must be a model for `Kernel`.

\tparam ET NT must be models for `RingNumberType`. Their default is  `K::RT`.


\cgalModels{PolytopeDistanceDTraits}

\sa `CGAL::Polytope_distance_d<Traits>`
\sa `CGAL::Polytope_distance_d_traits_2<K,ET,NT>`
\sa `CGAL::Polytope_distance_d_traits_3<K,ET,NT>`
\sa `PolytopeDistanceDTraits`

*/
template< typename K, typename ET, typename NT >
class Polytope_distance_d_traits_d {
public:

/// \name Types
/// @{

/*!
typedef to `K::Point_d`.
*/
typedef unspecified_type Point_d;

/*!
typedef to `K::Rep_tag`.
*/
typedef unspecified_type Rep_tag;

/*!
typedef to `K::RT`.
*/
typedef unspecified_type RT;

/*!
typedef to `K::FT`.
*/
typedef unspecified_type FT;

/*!
typedef to `K::Access_dimension_d`.
*/
typedef unspecified_type Access_dimension_d;

/*!
typedef to `K::Access_coordinates_begin_d`.
*/
typedef unspecified_type Access_coordinates_begin_d;

/*!
typedef to `K::Construct_point_d`.
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
Polytope_distance_d_traits_d( );

/*!
copy constructor.
*/
Polytope_distance_d_traits_d(
const Polytope_distance_d_traits_d<K,ET,NT>&);

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

}; /* end Polytope_distance_d_traits_d */
} /* end namespace CGAL */

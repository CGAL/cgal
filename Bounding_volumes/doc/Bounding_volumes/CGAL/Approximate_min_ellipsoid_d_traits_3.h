
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

The class `Approximate_min_ellipsoid_d_traits_3` is a traits class for
`CGAL::Approximate_min_ellipsoid_d<Traits>` using the
3-dimensional \cgal kernel. In order to use this class,
an exact number-type `ET` has to be provided which
`Approximate_min_ellipsoid_d<Traits>` will use for its internal
exact computations.

\tparam K must be a model for  `Kernel`.

\tparam ET must be a model for the concept
`EuclideanRing` with exact arithmetic operations, i.e., the type
\link AlgebraicStructureTraits::Is_exact `Algebraic_structure_traits<ET>::Is_exact` \endlink must be
`Tag_true`  (Examples of such a
number-type are `MP_Float`, `CORE::Expr`, and `Gmpq`.)

\cgalModels `ApproximateMinEllipsoid_d_Traits_d`

\sa `CGAL::Approximate_min_ellipsoid_d_traits_2<K,ET>`
\sa `CGAL::Approximate_min_ellipsoid_d_traits_d<K,ET>`
\sa `ApproximateMinEllipsoid_d_Traits_d`

*/
template< typename K, typename ET >
struct Approximate_min_ellipsoid_d_traits_3 {

/// \name Types
/// @{

/*!
`typedef double FT`. The kernel's number type
`K::RT` must be convertible to `double`.
*/
typedef unspecified_type FT;

/*!
typedef to the second template argument, `ET`.
*/
typedef unspecified_type ET;

/*!
`typedef K::Point_3 Point`
*/
typedef unspecified_type Point;

/*!
`typedef K::Cartesian_const_iterator_3 Cartesian_const_iterator`
*/
typedef unspecified_type Cartesian_const_iterator;

/// @}

}; /* end Approximate_min_ellipsoid_d_traits_3 */
} /* end namespace CGAL */

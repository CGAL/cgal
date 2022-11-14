
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

The class `Approximate_min_ellipsoid_d_traits_d` is a traits class for
`CGAL::Approximate_min_ellipsoid_d<Traits>` using the
d-dimensional \cgal kernel. In order to use this class,
an exact number-type `ET` has to be provided which
`Approximate_min_ellipsoid_d<Traits>` will use for its internal
exact computations.

\tparam K must be a model for concept `Kernel`.
\tparam ET must be a model for concept
`EuclideanRing` with exact arithmetic operations, i.e., the type
\link AlgebraicStructureTraits::Is_exact `Algebraic_structure_traits<ET>::Is_exact` \endlink must be
`CGAL::Tag_true` (Examples of such a
number-type are `MP_Float`, `CORE::Expr`, and `Gmpq`.)

\cgalModels `ApproximateMinEllipsoid_d_Traits_d`

\sa `CGAL::Approximate_min_ellipsoid_d_traits_2<K,ET>`
\sa `CGAL::Approximate_min_ellipsoid_d_traits_3<K,ET>`
\sa `ApproximateMinEllipsoid_d_Traits_d`

*/
template< typename K, typename ET >
struct Approximate_min_ellipsoid_d_traits_d {

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
`typedef K::Point_d Point`
*/
typedef unspecified_type Point;

/*!
`typedef K::Cartesian_const_iterator Cartesian_const_iterator`
*/
typedef unspecified_type Cartesian_const_iterator;

/// @}

}; /* end Approximate_min_ellipsoid_d_traits_d */
} /* end namespace CGAL */

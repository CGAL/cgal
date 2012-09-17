
namespace CGAL {

/*!
\ingroup PkgBoundingVolumes

The class `Approximate_min_ellipsoid_d_traits_2` is a traits class for 
`CGAL::Approximate_min_ellipsoid_d<Traits>` using the 
\f$ 2\f$-di-men-sional \cgal kernel. In order to use this class, 
an exact number-type `ET` has to be provided which 
`Approximate_min_ellipsoid_d<Traits>` will use for its internal 
exact computations. 

\tparam K must be a model for concept `Kernel`. 

\tparam ET must be a model for concept `RingNumberType` with exact arithmetic operations, i.e., the type 
`CGAL::Number_type_traits<ET>::Has_exact_ring_operations` must be 
`CGAL::Tag_true`. In addition, `ET` must be able to exactly 
represent any finite `double` value. (Examples of such a 
number-type are `CGAL::MP_Float`, `CORE::Expr`, and `CGAL::Gmpq`.) 

\models ::ApproximateMinEllipsoid_d_Traits_d 

\sa `CGAL::Approximate_min_ellipsoid_d_traits_3<K,ET>` 
\sa `CGAL::Approximate_min_ellipsoid_d_traits_d<K,ET>` 
\sa `ApproximateMinEllipsoid_d_Traits_d` 

*/
template< typename K, typename ET >
class Approximate_min_ellipsoid_d_traits_2 {
public:

/// \name Types 
/// @{

/*! 
`typedef double FT`. The kernel's number type 
`K::RT` must be convertible to `double`. 
*/ 
typedef Hidden_type FT; 

/*! 
typedef to the second template argument, `ET`. 
*/ 
typedef Hidden_type ET; 

/*! 
`typedef K::Point_2 Point` 
*/ 
typedef Hidden_type Point; 

/*! 
`typedef K::Cartesian_const_iterator_2 Cartesian_const_iterator` 
*/ 
typedef Hidden_type Cartesian_const_iterator; 

/// @}

}; /* end Approximate_min_ellipsoid_d_traits_2 */
} /* end namespace CGAL */

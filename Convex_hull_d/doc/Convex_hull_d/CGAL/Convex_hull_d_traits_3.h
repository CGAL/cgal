namespace CGAL {

/*!
\ingroup PkgConvexHullDRef

\deprecated This package is deprecated since the version 4.6 of \cgal. The package \ref PkgTriangulations should be used instead.

<A NAME="Index_anchor_0"></A>

The class `Convex_hull_d_traits_3` serves as a traits class for the class
`Convex_hull_d`. This is a traits class that adapts any
low-dimensional standard kernel model, e.g. `Homogeneous<RT>` or
`Cartesian<FT>` for the fixed 3-dimensional usage of
`Convex_hull_d`.

\cgalModels{ConvexHullTraits_d}

*/
template< typename R >
struct Convex_hull_d_traits_3 {

/// \name Creation
/// @{

/*!
default constructor.
*/
Convex_hull_d_traits_3();

/// @}

}; /* end Convex_hull_d_traits_3 */
} /* end namespace CGAL */

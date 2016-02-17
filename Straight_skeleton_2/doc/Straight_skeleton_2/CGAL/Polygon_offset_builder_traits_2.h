namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Classes

The class `Polygon_offset_builder_traits_2` provides a model for the 
`PolygonOffsetBuilderTraits_2` concept which is the traits class 
required by the `Polygon_offset_builder_2` algorithm class. The class 
\tparam Kernel a 2D \cgal Kernel such as the `Exact_predicates_inexact_constructions_kernel` (recommended)

It is unspecified which subset of the kernel is used in this traits class (and by extension in the builder class). This is to avoid restricting the choices in the implementation. 

\cgalModels `PolygonOffsetBuilderTraits_2`
\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`

\sa `Polygon_offset_builder_2`

*/
template< typename Kernel >
class Polygon_offset_builder_traits_2 {
public:

/// @}

}; /* end Polygon_offset_builder_traits_2 */
} /* end namespace CGAL */

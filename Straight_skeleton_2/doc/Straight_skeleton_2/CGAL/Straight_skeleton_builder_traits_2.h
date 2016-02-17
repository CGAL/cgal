namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Classes

The class `Straight_skeleton_builder_traits_2` provides a model for the 
`StraightSkeletonBuilderTraits_2` concept which is the traits class 
required by the `Straight_skeleton_builder_2` algorithm class.
\tparam Kernel a 2D \cgal Kernel,  such as the `Exact_predicates_inexact_constructions_kernel` (recommanded)

It is unspecified which subset of the kernel is used into the output sequence and the returned iterator will be equal to `out`. 

For any given input polygon, it in this traits class (and by extension in the builder class). This is to avoid restricting the choices in the implementation. 

\cgalModels `StraightSkeletonBuilderTraits_2`
\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`

\sa `CGAL::Straight_skeleton_builder_2<Gt,Ssds>` 

*/
template< typename Kernel >
class Straight_skeleton_builder_traits_2 {
public:

/// @}

}; /* end Straight_skeleton_builder_traits_2 */
} /* end namespace CGAL */

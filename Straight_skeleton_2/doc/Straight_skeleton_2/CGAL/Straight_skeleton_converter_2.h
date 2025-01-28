namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a new straight skeleton data structure with the same combinatorial and geometric data as `s`
using the items converter `ic` to convert the geometric embedding to the types of the target skeleton.

\tparam Target_skeleton must be a model of `StraightSkeleton_2`
\tparam Source_skeleton must be a model of `StraightSkeleton_2`
\tparam Items_converter must be a model of `StraightSkeletonItemsConverter_2`

\sa `CGAL::Straight_skeleton_items_converter_2<SrcSs,TgtSs,NTCV>`
\sa `CGAL::Straight_skeleton_converter_2<SrcSs,TgtSs,ItemsCV>`
*/
template <class Target_skeleton, class Source_skeleton, class Items_converter>
std::shared_ptr<Target_skeleton>
convert_straight_skeleton_2( Source_skeleton const& s,
                             Items_converted const& ic = Items_converter() );

/*!
\ingroup PkgStraightSkeleton2Auxiliary

The class `Straight_skeleton_converter_2` converts a straight skeleton instantiated using certain traits
into another straight skeleton instantiated using a different traits.

\tparam Source_skeleton_ type of the source straight skeleton
\tparam Target_skeleton_ type of the target straight skeleton
\tparam Items_converter_ a model `StraightSkeletonItemsConverter_2`. The default value of this parameter
                         is `Straight_skeleton_items_converter_2<Source_skeleton_,Target_skeleton_>`.

This conversion can be used to produce a straight skeleton using a kernel without exact constructions
(such as `Exact_predicates_inexact_constructions_kernel`) but input that skeleton
into `Polygon_offset_builder<Ss,Gt,Container>` instantiated
with a slower kernel (such as `Exact_predicates_exact_constructions_kernel`)
thus obtaining only simple offset polygons without paying the runtime overhead
of exact constructions for the straight skeleton itself.

\sa `convert_straight_skeleton_2()`
*/
template< typename Source_skeleton_, typename Target_skeleton_, typename Items_converter_ >
struct Straight_skeleton_converter_2 {

/// \name Types
/// @{

/*!
The `Source_skeleton_` template parameter corresponding to the source straight skeleton
*/
typedef Source_skeleton_ Source_skeleton;

/*!
The `Target_skeleton_` template parameter corresponding to the target straight skeleton
*/
typedef Target_skeleton_ Target_skeleton;

/*!
The `ItemsCvt` template parameter corresponding to the items converter
*/
typedef Items_converter_ Items_converter;

/// @}

/// \name Creation
/// @{

/*!
%Default constructor .
*/
Straight_skeleton_converter_2( const Items_converter& c = Items_converter() );

/// @}

/// \name Operations
/// @{

/*!
returns a new straight skeleton data structure with the same combinatorial and geometric data as `s` using the items converter to convert the geometric embeeding to the types of the target traits.
*/
std::shared_ptr<Target_skeleton> operator()( const Source_skeleton& s) const;

/// @}

}; /* end Straight_skeleton_converter_2 */


/*!
\ingroup PkgStraightSkeleton2Auxiliary

\cgalModels{StraightSkeletonItemsConverter_2}

`Straight_skeleton_items_converter_2` is a model of the `StraightSkeletonItemsConverter_2` concept

\tparam Source_skeleton_ type of the source straight skeleton
\tparam Target_skeleton_ type of the target straight skeleton
\tparam NT_converter_ a function object that must
                      provide `Target_skeleton_::Traits::FT operator()(Target_skeleton_::Traits::FT n)`
                      that converts `n` to an `Target_skeleton_::Traits::FT` which has the same value.
                      The default value of this parameter is `NT_converter<Source_skeleton_::Traits::FT, Target_skeleton_::Traits::FT>`.

\sa `convert_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_converter_2<SrcSs,TgtSs,ItemCV>`
*/
template< typename Source_skeleton_, typename Target_skeleton_, typename NT_converter_ >
struct Straight_skeleton_items_converter_2 {

/// \name Creation
/// @{

/*!
%Default constructor.
*/
Straight_skeleton_items_converter_2();

/// @}

}; /* end Straight_skeleton_items_converter_2 */

} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Functions

Converts a straight skeleton instantiated using certain traits into another straight skeleton instantiated using a different traits. 
returns a new straight skeleton data structure with the same combinatorial and geometric data as `s` using the items converter `ic` to convert the geometric embedding to the types of the target skeleton.

\tparam Target_skeleton type of the target straight skeleton
\tparam Source_skeleton type of the source straight skeleton
\tparam Items_converter a model of `StraightSkeletonItemsConverter_2`

\sa `StraightSkeletonItemsConverter_2`
\sa `Straight_skeleton_items_converter_2`
\sa `Straight_skeleton_converter_2`

*/
template <class Target_skeleton, class Source_skeleton, class Items_converter>
boost::shared_ptr<Target_skeleton, Source_skeleton, Items_converter>
convert_straight_skeleton_2( Source_skeleton const& s
, Items_converted const& ic = Items_converter()
);

/*!
\ingroup PkgStraightSkeleton2Auxiliary

The class `Straight_skeleton_converter_2` converts a straight skeleton instantiated using certain traits into another straight skeleton instantiated using a different traits. 


The third template parameter `ItemsConverter` is a function object that must 


\tparam SrcSs type of the source straight skeleton
\tparam TgtSs type of the target straight skeleton
\tparam ItemsCvt a model `StraightSkeletonItemsConverter_2`. The default value of this parameter is `Straight_skeleton_items_converter_2<SrcSs,TgtSs>`. 


This conversion can be used to produce a straight skeleton with the fast 
`Exact_predicates_inexact_constructions_kernel` but then input that skeleton 
into `Polygon_offset_builder<Ss,Gt,Container>` instantiated with the slower `Exact_predicates_exact_constructions_kernel` obtaining only simple offset 
polygons without paying the runtime overhead of exact constructions for the straight 
skeleton itself. 

\sa `StraightSkeletonItemsConverter_2` 
\sa `Straight_skeleton_items_converter_2` 
\sa `convert_straight_skeleton_2()`

*/
template< typename SrcSs, typename TgtSs, typename ItemsCvt >
class Straight_skeleton_converter_2 {
public:

/// \name Types 
/// @{

/*!
The `SrcSs` template parameter corresponding to the source straight skeleton. 
*/ 
SrcSs Source_skeleton; 

/*!
The `TgtSs` template parameter corresponding to the target straight skeleton. 
*/ 
TgtSs Target_skeleton; 

/*!
The `ItemsCvt` template parameter corresponding to the items converter. 
*/ 
ItemsCvt Items_converter; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor .
*/ 
Straight_skeleton_converter_2( Items_converter const& c ); 

/// @} 

/// \name Operations 
/// @{

/*!
returns a new straight skeleton data structure with the same combinatorial and geometric data as `s` using the items converter to convert the geometric embeeding to the types of the target traits. 
*/ 
boost::shared_ptr<Target_skeleton> operator()( Source_skeleton const& s) const; 

/// @}

}; /* end Straight_skeleton_converter_2 */


/*!
\ingroup PkgStraightSkeleton2Auxiliary

`Straight_skeleton_items_converter_2` is a model of the `StraightSkeletonItemsConverter_2` concept 

\tparam SrcSs type of the source straight skeleton
\tparam TgtSs type of the target straight skeleton
\tparam NTConverter  a function object that must 
provide `TgtSs:Traits::FT operator()(SrcSs::Traits::FT n)` that converts `n` to an 
`TgtSs::Traits::FT` which has the same value. The default value of this parameter is `NT_converter<SrcSs::Traits::FT, TgtSs::Traits::FT>`. 

\cgalModels `StraightSkeletonItemsConverter_2`
\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`
\cgalModels `Assignable`

\sa `Straight_skeleton_2`
\sa `Straight_skeleton_converter_2`
\sa `Straight_skeleton_builder_2`

*/
template< typename SrcSs, typename TgtSs, typename NTConverter >
class Straight_skeleton_items_converter_2 {
public:

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Straight_skeleton_items_converter_2<>(); 

/// @}

}; /* end Straight_skeleton_items_converter_2 */
} /* end namespace CGAL */


namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2

The class `Segment_Delaunay_graph_filtered_traits_2` provides a model for the 
`SegmentDelaunayGraphTraits_2` concept. 

The class `Segment_Delaunay_graph_filtered_traits_2` uses the filtering technique \cgalCite{cgal:bbp-iayed-01} 
to achieve traits for the `Segment_Delaunay_graph_2<Gt,DS>` 
class with efficient and exact predicates given an exact 
kernel `EK` and a filtering kernel `FK`. The geometric 
constructions associated provided by this class are equivalent 
to those provided by the traits class 
`Segment_Delaunay_graph_traits_2<CK,CM>`, which means that 
they may be inexact depending on the choice of the `CK` kernel. 

This class has six template parameters. 
\tparam CK is the construction kernel and it is the kernel 
that will be used for constructions.
\tparam CM must be `Field_with_sqrt_tag` or `Field_tag`
\tparam FK is the  filtering kernel; this kernel will be used for performing the 
arithmetic filtering for the predicates involved in the computation of 
the segment Delaunay graph. 
\tparam FM must be `Field_with_sqrt_tag` or `Field_tag`
\tparam EK is the exact kernel; this kernel will be used for computing the predicates if 
the filtering kernel fails to produce an answer. 
\tparam EM must be `Field_with_sqrt_tag` or `Field_tag`

The first, third and fifth template parameters must be models of the `Kernel` concept.
The second, fourth and sixth template parameters correspond to how 
predicates are evaluated. There are two predefined possible values for 
these parameters, namely `Field_with_sqrt_tag` and 
`Field_tag`. The first one must be used when the number type 
used in the representation supports the exact evaluation of signs of 
expressions involving all four basic operations and square roots, 
whereas the second one requires that only field operations are 
exact. Finally, in order to get exact constructions `CM` 
must be set to `Field_with_sqrt_tag` and the number type in 
`CK` must support operations involving divisions and square roots 
(as well as the other three basic operations of course). 
The way the predicates are evaluated is discussed in 
\cgalCite{b-ecvdl-96} and \cgalCite{cgal:k-reisv-04} (the geometric filtering 
part). 

The default values for the template parameters are as follows:
<ul>
 <li> `CM = Field_with_sqrt_tag` (it is assumed that `Cartesian<double>` or `Simple_cartesian<double>`  will be the entry for the template parameter `CK`), 
 <li> `EM = Field_tag`, 
 <li> `FK = Simple_cartesian<Interval_nt<false> >`, 
 <li> `FM = Field_with_sqrt_tag`. 
 <li> If the \sc{Gmp} package is  installed with \cgal, the template parameter `EK` has the default   value: `EK = Simple_cartesian<Gmpq>`, otherwise its   default value is   `EK = Simple_cartesian<Quotient<MP_Float> >`. 
</ul>

\cgalModels `SegmentDelaunayGraphTraits_2`
\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`
\cgalModels `Assignable`

\sa `Kernel` 
\sa `SegmentDelaunayGraphTraits_2` 
\sa `CGAL::Field_tag` 
\sa `CGAL::Field_with_sqrt_tag` 
\sa `CGAL::Segment_Delaunay_graph_2<Gt,DS>` 
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS>` 
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>` 

*/
template< typename CK, typename CM, typename EK, typename EM, typename FK, typename FM >
class Segment_Delaunay_graph_filtered_traits_2 {
public:

/// \name Types 
/// In addition to the types required by the
/// `SegmentDelaunayGraphTraits_2` concept the class
/// `Segment_Delaunay_graph_filtered_traits_2` defines the following
/// types:
/// @{

/*!

*/ 
typedef CGAL::Tag_true Intersections_tag; 

/*!

*/ 
typedef CK Kernel; 

/*!

*/ 
typedef CK Construction_kernel; 

/*!

*/ 
typedef FK Filtering_kernel; 

/*!

*/ 
typedef EK Exact_kernel; 

/*!

*/ 
typedef CM Method_tag; 

/*!

*/ 
typedef CM Construction_traits_method_tag; 

/*!

*/ 
typedef FM Filtering_traits_method_tag; 

/*!

*/ 
typedef EM Exact_traits_method_tag; 

/*!
A type for the segment Delaunay 
graph traits, where the kernel is `CK`. 
*/ 
typedef unspecified_type Construction_traits; 

/*!
A type for the segment Delaunay 
graph traits, where the kernel is `FK`. 
*/ 
typedef unspecified_type Filtering_traits; 

/*!
A type for the segment Delaunay 
graph traits, where the kernel is `EK`. 
*/ 
typedef unspecified_type Exact_traits; 

/// @}

}; /* end Segment_Delaunay_graph_filtered_traits_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2

The class `Segment_Delaunay_graph_filtered_traits_without_intersections_2` provides a model for the 
`SegmentDelaunayGraphTraits_2` concept. 

The class `Segment_Delaunay_graph_filtered_traits_without_intersections_2` uses the filtering technique \cgalCite{cgal:bbp-iayed-01} 
to achieve traits for the `Segment_Delaunay_graph_2<Gt,DS>` 
class with efficient and exact predicates given an exact 
kernel `EK` and a filtering kernel `FK`. The geometric 
constructions associated provided by this class are equivalent 
to those provided by the traits class 
`Segment_Delaunay_graph_traits_without_intersections_2<CK,CM>`, 
which means that they may be inexact, depending on the choice of the 
`CK` kernel. 

This class has six template parameters. 

\tparam CK is the construction kernel and it is the kernel 
that will be used for constructions.
\tparam CM must be `Field_with_sqrt_tag` or `Field_tag`
\tparam FK is the 
filtering kernel; this kernel will be used for performing the 
arithmetic filtering for the predicates involved in the computation of 
the segment Delaunay graph. 
\tparam FM must be `Field_with_sqrt_tag` or `Field_tag`
\tparam EK is the exact kernel; this kernel will be used for computing the predicates if 
the filtering kernel fails to produce an answer. 
\tparam EM must be `Field_with_sqrt_tag` or `Field_tag`

The first, third and fifth 
template parameters must be a models of the `Kernel` concept.
The second, fourth and sixth template parameters 
correspond to how predicates are evaluated. There are two predefined 
possible values for these parameters, namely `Field_with_sqrt_tag` 
and `Euclidean_ring_tag`. The first one must be used when the number 
type used in the representation supports the exact evaluation of signs 
of expressions involving all four basic operations and square roots, 
whereas the second requires the exact evaluation of signs of ring-type 
expressions, i.e., expressions involving only additions, subtractions 
and multiplications. Finally, in order to get exact constructions 
`CM` must be set to `Field_with_sqrt_tag` and the number type 
in `CK` must support operations involing divisions and square 
roots (as well as the other three basic operations of course). 
The way the predicates are evaluated is discussed in 
\cgalCite{b-ecvdl-96} and \cgalCite{cgal:k-reisv-04} (the geometric filtering 
part). 

The default values for the template parameters are as follows:
<ul>
<li> `CM = CGAL::Field_with_sqrt_tag` (it is assumed that 
`Cartesian<double>` or `Simple_cartesian<double>` 
will be the entry for the template parameter `CK`), 
<li> `EM = CGAL::Euclidean_ring_tag`, 
<li> `FK = CGAL::Simple_cartesian<CGAL::Interval_nt<false> >`, 
<li> `FM = CGAL::Field_with_sqrt_tag`. 
<li> If the \sc{Gmp} package is 
installed with \cgal, the template parameter `EK` has the default 
value: `EK = CGAL::Simple_cartesian<CGAL::Gmpq>`, otherwise its 
default value is 
`EK = CGAL::Simple_cartesian<CGAL::MP_Float>`. 
</ul>

\cgalModels `SegmentDelaunayGraphTraits_2`
\cgalModels `DefaultConstructible`
\cgalModels `CopyConstructible`
\cgalModels `Assignable`

\sa `Kernel` 
\sa `SegmentDelaunayGraphTraits_2` 
\sa `CGAL::Euclidean_ring_tag` 
\sa `CGAL::Field_with_sqrt_tag` 
\sa `CGAL::Segment_Delaunay_graph_2<Gt,DS>` 
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS>` 
\sa `CGAL::Segment_Delaunay_graph_traits_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>` 
\sa `CGAL::Segment_Delaunay_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>` 

*/
template< typename CK, typename CM, typename EK, typename EM, typename FK, typename FM >
class Segment_Delaunay_graph_filtered_traits_without_intersections_2 {
public:

/// \name Types 
/// In addition to the types required by the
/// `SegmentDelaunayGraphTraits_2` concept the class
/// `Segment_Delaunay_graph_filtered_traits_without_intersections_2`
/// defines the following types:
/// @{

/*!

*/ 
typedef CGAL::Tag_false Intersections_tag; 

/*!

*/ 
typedef CK Kernel; 

/*!

*/ 
typedef CK Construction_kernel; 

/*!

*/ 
typedef FK Filtering_kernel; 

/*!

*/ 
typedef EK Exact_kernel; 

/*!

*/ 
typedef CM Method_tag; 

/*!

*/ 
typedef CM Construction_traits_method_tag; 

/*!

*/ 
typedef FM Filtering_traits_method_tag; 

/*!

*/ 
typedef EM Exact_traits_method_tag; 

/*!
A type for the segment Delaunay 
graph traits, where the kernel is `CK`. 
*/ 
typedef unspecified_type Construction_traits; 

/*!
A type for the segment Delaunay 
graph traits, where the kernel is `FK`. 
*/ 
typedef unspecified_type Filtering_traits; 

/*!
A type for the segment Delaunay 
graph traits, where the kernel is `EK`. 
*/ 
typedef unspecified_type Exact_traits; 

/// @}

}; /* end Segment_Delaunay_graph_filtered_traits_without_intersections_2 */
} /* end namespace CGAL */

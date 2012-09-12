
namespace CGAL {

/*!
\ingroup PkgArrangement2

The traits class `Arr_polyline_traits_2` is a model of the `ArrangementTraits_2` 
concept. It handles piecewise linear curves, commonly referred to as 
polylines. Each polyline is a chain of segments, where each two neighboring 
segments in the chain share a common endpoint. The traits class exploits the 
functionality of the `SegmentTraits` template-parameter to handle the 
segments that comprise the polyline curves. 

The class instantiated for the template parameter `SegmentTraits` must 
be a model of the `ArrangementTraits_2` concept that handles line 
segments (e.g., `Arr_segment_traits_2<Kernel>` or 
`Arr_non_caching_segment_traits_2<Kernel>`, where the first 
alternative is recommended). 

The number type used by the injected segment traits should support exact 
rational arithmetic (that is, the number type should support 
the arithmetic operations \f$ +\f$, \f$ -\f$, \f$ \times\f$ and \f$ \div\f$ that should be 
carried out without loss of precision), in order to avoid robustness 
problems, although other inexact number types could be used at the user's 
own risk. 

\models ::ArrangementTraits_2 
\models ::ArrangementLandmarkTraits_2 
CONVERRORIsModel: CONVERROR 2 nested classes missing 

\sa `Arr_segment_traits_2<Kernel>` 
\sa `Arr_non_caching_segment_traits_2<Kernel>` 

*/
template< typename SegmentTraits >
class Arr_polyline_traits_2 {
public:

/// @}

}; /* end Arr_polyline_traits_2 */
} /* end namespace CGAL */

namespace CGAL {
namespace Box_intersection_d {
/*!
\ingroup PkgBoxIntersectionDClasses

This is the default traits class for the intersection algorithms for 
iso-oriented boxes. There are actually three versions depending on the 
type of `BoxHandle`; there is one if `BoxHandle` is a class type and there are 
two if `BoxHandle` is a pointer type, one for a mutable and one for a const 
pointer, respectively. 

This class implements the mapping from its `BoxHandle` argument to 
the `Box_parameter` type required in the 
`BoxIntersectionTraits_d` concept. In particular in the case where 
`BoxHandle` is a class type `B`, it defines 
`Box_parameter` to be of type `const B&`, while for the other 
cases it just uses the pointer type. 

<UL> 
<LI>`BoxHandle`: either a class type `B`, a pointer `B*`, or a 
const-pointer `const B*`, where `B` is a model of the 
`BoxIntersectionBox_d` concept. 
</UL> 

\cgalModels `BoxIntersectionTraits_d`

\sa \link PkgBoxIntersectionD_box_intersection_d `CGAL::box_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_d `CGAL::box_self_intersection_d()` \endlink
\sa \link PkgBoxIntersectionD_box_intersection_all_pairs_d `CGAL::box_intersection_all_pairs_d()` \endlink
\sa \link PkgBoxIntersectionD_box_self_intersection_all_pairs_d `CGAL::box_self_intersection_all_pairs_d()` \endlink
\sa `BoxIntersectionBox_d` 
\sa `CGAL::Box_intersection_d::Box_d<NT,int D,IdPolicy>` 

*/
template< typename BoxHandle >
class Box_traits_d {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Box_traits_d(); 

/// @}

}; /* end Box_traits_d */


/*!
\ingroup PkgBoxIntersectionDFunctions
*/
enum Setting  { COMPLETE, BIPARTITE };

/*!
\ingroup PkgBoxIntersectionDFunctions
*/
enum Topology { HALF_OPEN, CLOSED };


} /* Box_intersection_d */
} /* end namespace CGAL */

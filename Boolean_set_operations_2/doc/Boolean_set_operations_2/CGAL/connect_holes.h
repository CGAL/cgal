namespace CGAL {

/*!
\addtogroup boolean_connect_holes Functions on Polygon with Holes
\ingroup PkgBooleanSetOperations2
\anchor ref_bso_connect_holes 


*/

/// @{
/*!
Connects the holes of `pwh` with its outer boundary. This is done 
by locating the topmost vertex in each hole in the polygon with holes 
`pwh`, and connecting it by a vertical segment to the polygon 
feature located directly above it (a vertex or an edge of the outer 
boundary, or of another hole). The function produces an output 
sequence of points, which corresponds to the traversal of the vertices 
of the input polygon; this traversal starts from the outer boundary 
and moves to the holes using the auxiliary vertical segments that 
were added to connect the polygon with its holes. The value-type 
of `oi` is `Kernel::Point_2`. 
\pre The input polygon with holes `pwh` is bounded (namely it has a valid outer boundary). 
*/
template <class Kernel, class Container,
class OutputIterator>
OutputIterator
connect_holes(const Polygon_with_holes_2<Kernel,Container>& pwh,
              OutputIterator oi);
/// @}
} /* namespace CGAL */


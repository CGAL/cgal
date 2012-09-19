
/*!
\ingroup PkgKdsConcepts
\cgalconcept

This concept is for proxy objects which get notified when a kinetic regular triangulation changes. It inherits all the methods of `DelaunayTriangulationVisitor_3`. 

\models ::Kinetic::DelaunayTriangulationVisitor_3 

\hasModel `Kinetic::Regular_triangulation_visitor_base_3`
\hasModel `Kinetic::Regular_triangulation_event_log_visitor_3`

*/

class Kinetic::RegularTriangulationVisitor_3 {
public:

/// \name Operations 
/// @{

/*! 
The point defined by `Key` is about to move from the cell. 
*/ 
void pre_move(Key, Cell); 

/*! 
The point defined by `Key` just moved to the cell. 
*/ 
void post_move(Key, Cell); 

/// @}

}; /* end Kinetic::RegularTriangulationVisitor_3 */


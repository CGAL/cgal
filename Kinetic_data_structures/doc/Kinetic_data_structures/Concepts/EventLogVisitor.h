namespace Kinetic {

/*!
\ingroup PkgKdsConcepts
\cgalConcept

This concept is for visitors which maintain a text log of events. 

\cgalHasModel `CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_3`
\cgalHasModel `CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_2`
\cgalHasModel `CGAL::Kinetic::Regular_triangulation_event_log_visitor_3`
\cgalHasModel `CGAL::Kinetic::Sort_event_log_visitor` 
*/
class EventLogVisitor {
public:

/// \name Types 
/// @{

/*!
An iterator through strings defining the events that occurred. Each event is represented by a `std::string`. 
*/ 
typedef unspecified_type Event_iterator; 

/// @} 

/// \name Operations 
/// @{

/*!
Begin iterating through the events. 
*/ 
Event_iterator events_begin(); 

/*!

*/ 
Event_iterator events_end(); 

/// @}

}; /* end EventLogVisitor */

} /* end namespace Kinetic */

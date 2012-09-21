namespace Kinetic {

/*!
\ingroup PkgKdsConcepts
\cgalconcept

This concept is for visitors which maintain a text log of events. 

\hasModel `CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_3`
\hasModel `CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_2`
\hasModel `CGAL::Kinetic::Regular_trianglation_event_log_visitor_3`
\hasModel `CGAL::Kinetic::Sort_event_log_visitor` 
*/
class EventLogVisitor {
public:

/// \name Types 
/// @{

/*! 
An iterator through strings defining the events that occurred. Each event is represented by a `std::string`. 
*/ 
typedef Hidden_type Event_iterator; 

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

<<<<<<< .mine

} /* end namespace Kinetic */
=======
} /* end namespae KineticConcepts */
>>>>>>> .r72297

/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

The concept `StraightSkeletonBuilder_2_Visitor` describes the requirements of the visitor class required by the algorithm class `Straight_skeleton_builder_2<Gt,Ss,Visitor>` in its third template parameter. 

\cgalHasModel CGAL::Dummy_straight_skeleton_builder_2_visitor

\sa `CGAL::Straight_skeleton_builder_2<Gt,Ss>` 

*/

class StraightSkeletonBuilder_2_Visitor {
public:

/// \name Types 
/// @{

/*!
A constant handle to a straight skeleton halfedge. 
*/ 
typedef unspecified_type Halfedge_const_handle; 

/*!
A constant handle to a straight skeleton vertex. 
*/ 
typedef unspecified_type Vertex_const_handle; 

/// @} 

/// \name Operations 
/// @{

/*!
Called for each contour halfedge added to the skeleton). 
*/ 
void on_contour_edge_entered( Halfedge_const_handle const& ) const; 

/*!
Called before the initialization stage (when initial events are discovered) is started. 
*/ 
void on_initialization_started( std::size_t number_of_vertices ) const; 

/*!
Called after the events for contour vertex `v` have been discovered. `is_reflex` is true 
if this vertex has an internal angle \f$ > \pi\f$, or `is_degenerate` is true 
if the internal angle is \f$ = \pi\f$. 
*/ 
void on_contour_vertex_processed( Vertex_const_handle const& v 
, bool is_reflex 
, bool is_degenerate 
) const; 

/*!
Called after an edge event for nodes `node0` and `node1` has been discovered 
and put on the queue for later processing. 
*/ 
void on_edge_event_created( Vertex_const_handle const& node0 
, Vertex_const_handle const& node1 ) const ; 

/*!
Called after a slipt event for node `node` has been discovered 
and put on the queue for later processing. 
*/ 
void on_split_event_created( Vertex_const_handle const& node ) const ; 

/*!
Called after a pseudo slipt event for nodes `node0` and `node1` has been discovered 
and put on the queue for later processing. 
*/ 
void on_pseudo_split_event_created( Vertex_const_handle const& node0 
, Vertex_const_handle const& node1 ) const ; 

/*!
Called after all initial events have been discovered. 
*/ 
void on_initialization_finished() const; 

/*!
Called before the propagation stage (when events are poped off the queue and processed) 
is started. 
*/ 
void on_propagation_started() const; 

/*!
Called after an anhiliation event for nodes `node0` and `node1` has been processed. 
A new skeleton edge between these nodes has been added. 
*/ 
void on_anihiliation_event_processed ( Vertex_const_handle const& node0 
, Vertex_const_handle const& node1 
) const; 

/*!
Called after an edge for nodes `seed0` and `seed1` has been processed. 
Skeleton vertex `newnode` and edges from `node0` to `newnode` and `node1` to `newnode` 
has been added. 
*/ 
void on_edge_event_processed( Vertex_const_handle const& seed0 
, Vertex_const_handle const& seed1 
, Vertex_const_handle const& newnode 
) const; 

/*!
Called after a split event for node `seed` has been processed. 
Skeleton vertices `newnode0` and `newnode1` have been added. 
An skeleton edge from `seed` to `newnode0` has been added. 
In the final skeleton, `newnode1` is removed and only `newnode0` remains. 
*/ 
void on_split_event_processed( Vertex_const_handle const& seed 
, Vertex_const_handle const& newnode0 
, Vertex_const_handle const& newnode1 
) const; 

/*!
Called after a pseudo split event for nodes `seed0` and `seed1` has been processed. 
Skeleton vertices `newnode0` and `newnode1` have been added. 
Skeleton edges from `seed0` to `newnode0` and `seed1` to `newnode1` has been added. 
*/ 
void on_pseudo_split_event_processed( Vertex_const_handle const& seed0 
, Vertex_const_handle const& seed1 
, Vertex_const_handle const& newnode0 
, Vertex_const_handle const& newnode1 
) const; 

/*!
Called after vertex `v` has been marked as already processed. 
*/ 
void on_vertex_processed( Vertex_const_handle const& v ) const; 

/*!
Called after all events have been processed. 
*/ 
void on_propagation_finished() const; 

/*!
Called when the skeleton clean up (when multiple nodes are merged) is started. 
*/ 
void on_cleanup_started() const; 

/*!
Called when clean up finished. 
*/ 
void on_cleanup_finished() const; 

/*!
Called when the algorithm terminated. 
`finished_ok` is false if it terminated before completion or the resulting skeleton 
was found to be invalid. 
*/ 
void on_algorithm_finished ( bool finished_ok ) const; 

/*!
Called whenever an error was detected. 
`msg` is whatever error message accompanies the error. This pointer can be `null`. 
*/ 
void on_error( char const* msg ) const; 

/// @}

}; /* end StraightSkeletonBuilder_2_Visitor */

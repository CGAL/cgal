/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

The concept `StraightSkeletonBuilder_2_Visitor` describes the requirements of the visitor class
required by the algorithm class `CGAL::Straight_skeleton_builder_2` in its third template parameter.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Dummy_straight_skeleton_builder_2_visitor}
\cgalHasModelsEnd

\sa `CGAL::Straight_skeleton_builder_2`
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
void on_contour_edge_entered( const Halfedge_const_handle& ) const;

/*!
Called before the initialization stage (when initial events are discovered) is started.
*/
void on_initialization_started( std::size_t number_of_vertices ) const;

/*!
Called after the initial events involving the contour vertex `v` have been discovered.
*/
void on_initial_events_collected ( const Vertex_const_handle& v,
                                   bool is_reflex,
                                   bool is_degenerate) const;

/*!
Called after an edge event for nodes `node0` and `node1` has been discovered
and put on the queue for later processing.
*/
void on_edge_event_created( const Vertex_const_handle& node0,
                            const Vertex_const_handle& node1 ) const ;

/*!
Called after a slipt event for node `node` has been discovered
and put on the queue for later processing.
*/
void on_split_event_created( const Vertex_const_handle& node ) const ;

/*!
Called after a pseudo slipt event for nodes `node0` and `node1` has been discovered
and put on the queue for later processing.
*/
void on_pseudo_split_event_created( const Vertex_const_handle& node0,
                                    const Vertex_const_handle& node1 ) const ;

/*!
Called after all initial events have been discovered.
*/
void on_initialization_finished() const;

/*!
Called before the propagation stage (when events are popped off the queue and processed)
is started.
*/
void on_propagation_started() const;

/*!
Called after an annihilation event for nodes `node0` and `node1` has been processed.
A new skeleton edge between these nodes has been added.
*/
void on_anihiliation_event_processed ( const Vertex_const_handle& node0,
                                       const Vertex_const_handle& node1 ) const;

/*!
Called after an edge for nodes `seed0` and `seed1` has been processed.
Skeleton vertex `newnode` and edges from `node0` to `newnode` and `node1` to `newnode`
has been added.
*/
void on_edge_event_processed( const Vertex_const_handle& seed0,
                              const Vertex_const_handle& seed1,
                              const Vertex_const_handle& newnode ) const;

/*!
Called after a split event for node `seed` has been processed.
Skeleton vertices `newnode0` and `newnode1` have been added.
An skeleton edge from `seed` to `newnode0` has been added.
In the final skeleton, `newnode1` is removed and only `newnode0` remains.
*/
void on_split_event_processed( const Vertex_const_handle& seed,
                               const Vertex_const_handle& newnode0,
                               const Vertex_const_handle& newnode1 ) const;

/*!
Called after a pseudo split event for nodes `seed0` and `seed1` has been processed.
Skeleton vertices `newnode0` and `newnode1` have been added.
Skeleton edges from `seed0` to `newnode0` and `seed1` to `newnode1` has been added.
*/
void on_pseudo_split_event_processed( const Vertex_const_handle& seed0,
                                      const Vertex_const_handle& seed1,
                                      const Vertex_const_handle& newnode0,
                                      const Vertex_const_handle& newnode1) const;

/*!
Called after vertex `v` has been marked as already processed.
*/
void on_vertex_processed( const Vertex_const_handle& v ) const;

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

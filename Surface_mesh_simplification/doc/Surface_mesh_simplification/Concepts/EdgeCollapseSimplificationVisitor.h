
/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalconcept

The concept `EdgeCollapseSimplificationVisitor` describes the requirements for the <I>visitor object</I> which is used to track the edge collapse simplification algorithm. 

The several callbacks given as member functions in the visitor are called from certain places in the algorithm implementation. 

*/

class EdgeCollapseSimplificationVisitor {
public:

/// \name Types 
/// @{

/*! 
The type of the surface to simplify. Must be a model of the `EdgeCollapsableMesh` concept. 
*/ 
typedef Hidden_type ECM; 

/*! 
A field type representing the collapse cost 
*/ 
typedef Hidden_type FT; 

/*! 
The type of the edge profile cache. Must be a model of the `EdgeProfile` concept. 
*/ 
typedef Hidden_type Profile; 

/*! 
The point type for the surface vertex. Must be a model of `Point_3`. 
*/ 
typename CGAL::halfedge_graph_traits<ECM>::Point Point; 

/*! 
An integer type representing the number of edges 
*/ 
typedef Hidden_type size_type; 

/// @} 

/// \name Operations 
/// @{

/*! 
Called before the algorithm starts. 
*/ 
void OnStarted( ECM& surface ); 

/*! 
Called after the algorithm finishes. 
*/ 
void OnFinished ( ECM& surface ) ; 

/*! 
Called when the `StopPredicate` returned `true` 
(but not if the algorithm terminates because the surface could not be simplified any further). 

*/ 
void OnStopConditionReached( ECM& surface ) ; 

/*! 
Called during the <I>collecting phase</I> (when a cost is assigned to the edges), 
for each edge collected. 

*/ 
void OnCollected( Profile const& profile, boost::optional<FT> cost ); 

/*! 
Called during the <I>processing phase</I> (when edges are collapsed), 
for each edge that is selected. 

This method is called before the algorithm checks 
if the edge is collapsable. 

`cost` indicates the current collapse cost for the edge. 
If absent (meaning that it could not be computed) 
the edge will not be collapsed. 

`initial_count` and `current_count` refer to 
the number of edges. 

*/ 
void OnSelected( Profile const& profile 
, boost::optional<FT> cost 
, size_type initial_count 
, size_type current_count 
); 

/*! 
Called when an edge is about to be collapsed and replaced by a vertex 
whose position is `*placement`. 

If `placement` is absent (meaning that it could not be computed) 
the edge will not be collapsed. 

*/ 
void OnCollapsing( Profile const& profile 
, boost::optional<Point> placement 
); 

/*! 
Called for each selected edge which cannot be 
collapsed because doing so would change the topological 
type of the surface (turn it into a non-manifold 
for instance). 

*/ 
void OnNonCollapsable( Profile const& profile ); 

/// @}

}; /* end EdgeCollapseSimplificationVisitor */


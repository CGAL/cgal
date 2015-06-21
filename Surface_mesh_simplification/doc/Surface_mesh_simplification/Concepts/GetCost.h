
/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `GetCost` describes the requirements for the <I>policy function object</I> 
which gets the <I>collapse cost</I> of an edge. 

The cost returned is a `boost::optional` value (i.e.\ it can be absent). 
An <I>absent</I> cost indicates that the edge should not be collapsed. 
This could be the result of a computational limitation (such as overflow), 
or can be intentionally returned to prevent the edge from being collapsed. 

\cgalRefines `DefaultConstructible` 
\cgalRefines `CopyConstructible` 

\cgalHasModel `CGAL::Surface_mesh_simplification::Edge_length_cost<ECM>`
\cgalHasModel `CGAL::Surface_mesh_simplification::LindstromTurk_cost<ECM>`

*/

class GetCost {
public:


/// \name Operations 
/// @{

/*!
Computes and returns the cost of collapsing the edge (represented by its profile), 
using the calculated placement. 
\tparam Profile must be a model of `EdgeProfile`.
*/ 
  template <class Profile>
  boost::optional<typename Profile::FT> operator()( Profile const& edge_profile 
                                                    , boost::optional<typename Profile::Point> const& placement ) const; 

/// @}

}; /* end GetCost */


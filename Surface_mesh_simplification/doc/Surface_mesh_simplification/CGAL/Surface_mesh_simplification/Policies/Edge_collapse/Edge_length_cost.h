
namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Edge_length_cost` is a model for the `GetCost` concept,
which computes the collapse cost as the squared length of the edge. 

\tparam ECM is the type of surface being simplified, and must be a model of the `EdgeCollapsableMesh` concept. 


\models ::GetCost 

*/
template< typename ECM >
class Edge_length_cost {
public:

/// \name Creation 
/// @{

/*! 
Default constructor 
*/ 
Edge_length_cost<ECM>(); 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns the <I>collapse cost</I> as the squared distance between the points 
of the source and target vertices (that is, `profile.p0()` and `profile.p1()`. 

The `placement` argument is ignored. 

*/ 
result_type operator()( Profile const& profile 
, boost::optional<Point> const& placement ) const; 

/// @}

}; /* end Surface_mesh_simplification::Edge_length_cost */
} /* namespace Surface_mesh_simplification */
} /* end namespace CGAL */

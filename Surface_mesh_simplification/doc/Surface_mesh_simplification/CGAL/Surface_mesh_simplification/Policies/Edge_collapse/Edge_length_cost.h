
namespace CGAL {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Surface_mesh_simplification::Edge_length_cost` provides a model for the `GetCost` concept. 
It has one template argument: the type of surface being simplified, 
which must be a model of the `EdgeCollapsableMesh` concept. 
It computes the collapse cost as the squared length of the edge. 

\models ::GetCost 

*/
template< typename ECM >
class Surface_mesh_simplification::Edge_length_cost {
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
} /* end namespace CGAL */

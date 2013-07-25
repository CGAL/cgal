namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The function object class `Mesh_edge_criteria_3` is a model of `MeshEdgeCriteria_3`. It 
provides a bound for the size criterion. 

\cgalModels `MeshEdgeCriteria_3`

\sa `MeshCriteria_3` 
\sa `CGAL::Mesh_criteria_3<Tr>` 
\sa `MeshDomainField_3` 

*/
template< typename Tr >
class Mesh_edge_criteria_3 {
public:

/// \name Types 
/// @{

/*!
Numerical type. 
*/ 
typedef Tr::Geom_traits::FT FT; 

/// @} 

/// \name Creation 
/// @{

/*!
Returns an object to serve as criteria for edges. 
The argument `length_bound` is an upper bound 
for the length of the edges which are used to discretize the curve segments. 
Note that if one parameter is set to 0, then its corresponding criteria is ignored. 
*/ 
Mesh_edge_criteria_3( 
FT length_bound); 

/*!
Returns an object to serve as criteria for edges. The type `SizingField` 
must be a model of concept `MeshDomainField_3`. The behavior and semantic of the argument are the same 
as above, except that the length 
parameter is a functional instead of a constant. 
*/ 
template <class SizingField> 
Mesh_edge_criteria_3( 
SizingField length_bound); 

/// @}

}; /* end Mesh_edge_criteria_3 */
} /* end namespace CGAL */

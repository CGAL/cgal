namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The class `Mesh_cell_criteria_3` is a model of `MeshCellCriteria_3`. It provides, 
for the mesh tetrahedra, 
a uniform shape criteria 
and a sizing field which may be a uniform or variable field. 

\tparam Tr must be identical to the nested type 
`Triangulation` of the instance used as model of 
`MeshComplex_3InTriangulation_3`. 

\cgalModels `MeshCellCriteria_3` 

\sa `MeshCriteria_3` 
\sa `CGAL::Mesh_criteria_3<Tr>` 
\sa `CGAL::make_mesh_3()` 

*/
template< typename Tr >
class Mesh_cell_criteria_3 {
public:

/// \name Types 
/// @{

/*!
Numerical type 
*/ 
typedef Tr::FT FT; 

/// @} 

/// \name Creation 
/// @{

/*!
Returns an object to serve as default criteria for cells. The argument 
`radius_edge_bound` is the upper bound for the radius-edge ratio 
of the tetrahedra. The argument `radius_bound` is a uniform upper bound 
for the circumradii of the tetrahedra in the mesh. See 
section \ref introsecparam for further details. 
Note that if one parameter is set to 0, then its corresponding criteria is ignored. 
*/ 
Mesh_cell_criteria_3( 
FT radius_edge_bound, 
FT radius_bound); 

/*!
Returns an object to serve as default criteria for facets. The type `SizingField` must 
be a model of the concept `MeshDomainField_3`. The behavior and semantic of the arguments are the same 
as above, except that the radius bound parameter is a functional instead of a constant. 
*/ 
template<class SizingField> 
Mesh_cell_criteria_3( 
FT radius_edge_bound, 
SizingField radius_bound); 

/// @}

}; /* end Mesh_cell_criteria_3 */
} /* end namespace CGAL */

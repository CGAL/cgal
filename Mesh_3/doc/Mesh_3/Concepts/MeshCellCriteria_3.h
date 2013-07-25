/*!
\ingroup PkgMesh_3Concepts
\cgalConcept

The Delaunay refinement process involved in the 
template functions `make_mesh_3()` and `refine_mesh_3()` 
is guided by a set of elementary refinement criteria 
that concern either mesh tetrahedra or surface facets. 
The concept `MeshCellCriteria_3` describes the types that 
handle the refinement criteria for mesh tetrahedra. 

\cgalHasModel `CGAL::Mesh_cell_criteria_3<Tr>`

\sa `MeshEdgeCriteria_3` 
\sa `MeshFacetCriteria_3` 
\sa `MeshCriteria_3` 
\sa `CGAL::make_mesh_3()` 
\sa `CGAL::refine_mesh_3()` 

*/

class MeshCellCriteria_3 {
public:

/// \name Types 
/// @{

/*!
Handle type for the cells of the 
triangulation. Must match the `Cell_handle` type in the 
triangulation type used by the mesh generation function. 
*/ 
typedef unspecified_type Cell_handle; 

/*!
Type representing the quality of a 
cell. Must be a model of CopyConstructible and 
LessThanComparable. Between two cells, the one which has the lower 
quality must have the lower `Cell_quality`. 
*/ 
typedef unspecified_type Cell_quality; 

/*!
Type representing if a cell is bad or not. Must 
be convertible to `bool`. If it converts to `true` then the cell is bad, otherwise 
the cell is good with regard to the criteria. 

In addition, an object of this type must contain an object of type 
`Cell_quality` if it represents 
a bad cell. `Cell_quality` must be accessible by `operator*()`. 
Note that `boost::optional<Cell_quality>` is a natural model of this concept. 
*/ 
typedef unspecified_type Is_cell_bad; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns `Is_cell_bad` value of cell `c`. 
*/ 
Is_cell_bad operator()(Cell_handle c); 

/// @}

}; /* end MeshCellCriteria_3 */

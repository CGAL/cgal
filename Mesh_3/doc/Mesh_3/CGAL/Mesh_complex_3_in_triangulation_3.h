namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The class `Mesh_complex_3_in_triangulation_3` implements a data structure 
to store the 3D restricted Delaunay triangulation used by a mesh 
generation process. 

This class is a model of the concept 
`MeshComplexWithFeatures_3InTriangulation_3`. 


\tparam Tr can be instantiated with any 3D 
regular triangulation of \cgal provided that its 
vertex and cell base class are models of the concepts 
`MeshVertexBase_3` and `MeshCellBase_3`, respectively. 

\tparam  CornerIndex is the type of the indices for corners. It must match the `Corner_index` of the model 
of the `MeshDomainWithFeatures_3` concept used for mesh generation. 

\tparam CurveSegmentIndex is the type of the indices for curves segments. 
It must match the `Curve_segment_index` types of the model 
of the `MeshDomainWithFeatures_3` concept used for mesh generation. 

Those two last template parameters defaults to `int`, so that they can be ignored 
if the domain used for mesh generation does not include 0 and 1-dimensionnal features (i.e 
is a model of the concept `MeshDomain_3`). 

\cgalModels `MeshComplexWithFeatures_3InTriangulation_3` 

\sa `CGAL::make_mesh_3()` 
\sa `CGAL::refine_mesh_3()` 
\sa `MeshComplex_3InTriangulation_3` 
\sa `MeshComplexWithFeatures_3InTriangulation_3` 
\sa `MeshCellBase_3`, 
\sa `MeshVertexBase_3` 

*/
template< typename Tr, typename CornerIndex, typename CurveSegmentIndex >
class Mesh_complex_3_in_triangulation_3 {
public:

/// \name Types 
/// @{

/*!
%Index type. 
*/ 
typedef Tr::Vertex::Index Index; 

/*!
Surface index type. 
*/ 
typedef Tr::Cell::Surface_patch_index Surface_patch_index; 

/*!
Subdomain index type. 
*/ 
typedef Tr::Cell::Subdomain_index Subdomain_index; 

/*!
 Corner index type. 
*/ 
typedef CornerIndex Corner_index; 

/*!
Curve segment index type. 
*/ 
typedef CurveSegmentIndex Curve_segment_index; 

/// @} 

/// \name Operations 
/// @{

/*!
Outputs the mesh to `os` 
in medit format. 
*/ 
void output_to_medit(std::ofstream& os); 

/**
 * Outputs the outer boundary of the entire domain with facets oriented outward.
 */
std::ostream& output_boundary_to_off(std::ostream& out) const;

/**
 * Outputs the outer boundary of the selected subdomain with facets oriented outward.
 */
std::ostream& output_boundary_to_off(std::ostream& out, Subdomain_index subdomain) const;

/**
 * Outputs the surface facets with a consistent orientation at the interface of two subdomains.
 */
std::ostream& output_facets_in_complex_to_off(std::ostream& out) const;

/// @}

}; /* end Mesh_complex_3_in_triangulation_3 */
} /* end namespace CGAL */

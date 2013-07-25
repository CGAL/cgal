
/*!
\ingroup PkgSurfaceSubdivisionMethods3Concepts
\cgalConcept

Required member functions for the `Sqrt3Mask_3` concept. This 
policy concept of geometric computations is used in 
`CGAL::Subdivision_method_3::Sqrt3<Polyhedron_3, Mask>`. 

\cgalHasModel `CGAL::Sqrt3_mask_3<Polyhedron_3>`

\sa `CGAL::Subdivision_method_3`

*/

class Sqrt3Mask_3 {
public:

/// \name Operations 
/// @{

/*!

computes the subdivided point `pt` based on the neighborhood 
of the facet `f`. 
*/ 
void facet_node(Facet_handle f, Point_3& pt); 

/*!

computes the subdivided point `pt` based on the neighborhood 
of the vertex `v`. 
*/ 
void vertex_node(Vertex_handle v, Point& pt); 

/// @}

}; /* end Sqrt3Mask_3 */


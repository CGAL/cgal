
/*!
\ingroup PkgSurfaceSubdivisionMethods3Concepts
\cgalConcept

Required member functions for the `Sqrt3Mask_3` concept. This 
policy concept of geometric computations is used in 
`CGAL::Subdivision_method_3::Sqrt3<PolygonMesh, Mask>`. 

\cgalRefines `SubdivisionMask_3`
\cgalHasModel `CGAL::Sqrt3_mask_3<PolygonMesh>`

\sa `CGAL::Subdivision_method_3`

*/

class Sqrt3Mask_3 {
public:

/// \name Operations 
/// @{

/*! Constructor
 */

Sqrt3Mask_3(TriangleMesh& m);

/*!
computes the subdivided point `pt` based on the neighborhood 
of the facet `f`. 
*/ 
void facet_node(face_descriptor f, Point_3& pt); 

/*!
computes the subdivided point `pt` based on the neighborhood 
of the vertex `v`. 
*/ 
void vertex_node(vertex_descriptor v, Point& pt); 

/// @}

}; /* end Sqrt3Mask_3 */


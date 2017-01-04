
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

/// \name Types 
/// @{

/*!  The polygon mesh must be triangulated.

*/ 
  typedef unspecified_type PolygonMesh;

/// @}


/// \name Operations 
/// @{

/*! Constructor
 */

Sqrt3Mask_3(PolygonMesh& pmesh);

/*!
computes the subdivided point `pt` based on the neighborhood 
of the face `fd`. 
*/ 
void facet_node(face_descriptor fd, Point_3& pt); 

/*!
computes the subdivided point `pt` based on the neighborhood 
of the vertex `vd`. 
*/ 
void vertex_node(vertex_descriptor vd, Point& pt); 

/// @}

}; /* end Sqrt3Mask_3 */


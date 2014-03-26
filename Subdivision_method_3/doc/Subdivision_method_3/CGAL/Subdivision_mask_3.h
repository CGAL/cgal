
namespace CGAL {

/*!
\ingroup PkgSurfaceSubdivisionMethods3

A stencil determines a source neighborhood 
whose points contribute to the position of a refined point. 
The geometry mask of a stencil specifies 
the computation on the nodes of the stencil. 
`CatmullClark_mask_3` implements the geometry masks of 
Catmull-Clark subdivision on a `CGAL::Polyhedron_3<Cartesian>`. 

\tparam Polyhedron_3 must be a `CGAL::Polyhedron_3`
instantiated  with a %Cartesian kernel, which defines the `Point_3` for the vertices. 

\image html CCBorderMask.png
\image latex CCBorderMask.png

\cgalModels `PQQMask_3`

\sa `CGAL::Subdivision_method_3`

*/
template< typename Polyhedron_3 >
class CatmullClark_mask_3 {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
CatmullClark_mask_3<Polyhedron_3>(); 

/// @} 

/// \name Stencil functions 
/// @{

/*!

computes the Catmull-Clark facet-point `pt` of the facet `f`. 

*/ 
void facet_node(Facet_handle f, Point_3& pt); 

/*!

computes the Catmull-Clark edge-point `pt` of the edge `e`. 

*/ 
void edge_node(Edge_handle e, Point_3& pt); 

/*!

computes the Catmull-Clark vertex-point `pt` of the vertex `v`. 

*/ 
void vertex_node(Vertex_handle v, Point_3& pt); 

/*!

computes the Catmull-Clark edge-point `ept` and the 
Catmull-Clark vertex-point `vpt` of the border edge `e`. 

*/ 
void border_node(Halfedge_handle e, Point_3& ept, Point_3& vpt); 

/// @}

}; /* end CatmullClark_mask_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSurfaceSubdivisionMethods3

A stencil determines a source neighborhood 
whose points contribute to the position of a refined point. 
The geometry mask of a stencil specifies 
the computation on the nodes of the stencil. 
`DooSabin_mask_3` implements the geometry masks of 
Doo-Sabin subdivision on a `Polyhedron_3<Cartesian>`. 

\tparam Polyhedron_3 must be a `CGAL::Polyhedron_3`
instantiated  with a %Cartesian kernel, which defines the `Point_3` for the vertices. 

\image html DSCornerMask.png
\image latex DSCornerMask.png

\cgalModels `DQQMask_3`

\sa `CGAL::Subdivision_method_3`

*/
template< typename Polyhedron_3 >
class DooSabin_mask_3 {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
DooSabin_mask_3<Polyhedron_3>(); 

/// @} 

/// \name Stencil functions 
/// @{

/*!

computes the Doo-Sabin point `pt` of the vertex pointed 
by the halfedge `he`. 

*/ 
void corner_node(Halfedge_handle he, Point_3& pt); 

/// @}

}; /* end DooSabin_mask_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSurfaceSubdivisionMethods3

A stencil determines a source neighborhood 
whose points contribute to the position of a refined point. 
The geometry mask of a stencil specifies 
the computation on the nodes of the stencil. 
`Loop_mask_3` implements the geometry masks of 
Loop subdivision on a triangulated `Polyhedron_3<Cartesian>`. 

\tparam Polyhedron_3 must be a `CGAL::Polyhedron_3`
instantiated  with a %Cartesian kernel, which defines the `Point_3` for the vertices. 

\image html LoopBorderMask.png
\image latex LoopBorderMask.png

\cgalModels `PTQMask_3`

\sa `CGAL::Subdivision_method_3`

*/
template< typename Polyhedron_3 >
class Loop_mask_3 {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Loop_mask_3<Polyhedron_3>(); 

/// @} 

/// \name Stencil functions 
/// @{

/*!

computes the Loop edge-point `pt` of the edge `e`. 

*/ 
void edge_node(Edge_handle e, Point_3& pt); 

/*!

computes the Loop vertex-point `pt` of the vertex `v`. 

*/ 
void vertex_node(Vertex_handle v, Point_3& pt); 

/*!

computes the Loop edge-point `ept` and the 
Loop vertex-point `vpt` of the border edge `e`. 

*/ 
void border_node(Halfedge_handle e, Point_3& ept, Point_3& vpt); 

/// @}

}; /* end Loop_mask_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSurfaceSubdivisionMethods3

A stencil determines a source neighborhood 
whose points contribute to the position of a refined point. 
The geometry mask of a stencil specifies 
the computation on the nodes of the stencil. 
`Sqrt3_mask_3` implements the geometry masks of 
\f$ \sqrt{3}\f$ subdivision on a triangulated 
`CGAL::Polyhedron_3<Cartesian>`. 

\tparam Polyhedron_3 must be a `CGAL::Polyhedron_3`
instantiated  with a %Cartesian kernel, which defines the `Point_3` for the vertices. 

\cgalModels `Sqrt3Mask_3`

\sa `CGAL::Subdivision_method_3`

*/
template< typename Polyhedron_3 >
class Sqrt3_mask_3 {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Sqrt3_mask_3<Polyhedron_3>(); 

/// @} 

/// \name Stencil functions 
/// @{

/*!

computes the \f$ \sqrt{3}\f$ facet-point `pt` of the facet `f`. 

*/ 
void facet_node(Facet_handle f, Point_3& pt); 

/*!

computes the \f$ \sqrt{3}\f$ vertex-point `pt` of the vertex `v`. 

*/ 
void vertex_node(Vertex_handle v, Point& pt); 

/// @}

}; /* end Sqrt3_mask_3 */
} /* end namespace CGAL */

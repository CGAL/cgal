
namespace CGAL {
/// The namespace containing the subdivision methods.
namespace Subdivision_method_3 {

/*!
\addtogroup PkgSurfaceSubdivisionMethods3Functions

A subdivision method recursively refines a coarse mesh and 
generates an ever closer approximation to a smooth surface. 
`Subdivision_method_3` consists of four subdivision methods 
and their refinement hosts. Each refinement host is a template 
function of a polyhedron class and a 
geometry policy class. It refines the connectivity of the 
control mesh and computes the geometry of the refined mesh. 
The geometry computation is dedicated to the custom 
geometry policy. A geometry policy consists of functions 
that compute the new point based on the subdivision stencil. 
A stencil defines the footprint (a submesh of the control mesh) 
of a new point. 

The four supported refinement hosts are the 
primal quadrilateral quadrisection (PQQ), 
the primal triangle quadrisection (PTQ), 
the dual quadrilateral quadrisection (DQQ), 
and the \f$ \sqrt{3}\f$ triangulation. 
These refinements are respectively used in 
Catmull-Clark, Loop, Doo-Sabin and \f$ \sqrt{3}\f$ subdivision. 

\cgalHeading{Refinement Host}

A refinement host is a template function of 
a polyhedron class and a geometry mask class. It refines 
the input polyhedron, and computes new points through 
the geometry masks. 
`Subdivision_method_3` supports four refinement hosts: 
`PQQ`, `PTQ`, `DQQ` and `Sqrt3`. 

\image html RefSchemes.png
\image latex RefSchemes.png

\cgalHeading{Example}

This example program subdivides a polyhedral mesh with 
Catmull-Clark subdivision. 

\cgalExample{Subdivision_method_3/CatmullClark_subdivision.cpp} 


\sa `CGAL::CatmullClark_mask_3<Polyhedron_3>`
\sa `CGAL::Loop_mask_3<Polyhedron_3>`
\sa `CGAL::Sqrt3_mask_3<Polyhedron_3>`
*/
/// @{

/*!

applies the PQQ refinement on the control mesh `p` `step` times. 
The geometry of the refined mesh is computed by the geometry policy `mask`. 
This function overwrites the control mesh `p` with the refined mesh. 
*/ 

template <class Polyhedron_3, template <typename> class Mask> 
void PQQ(Polyhedron_3& p, Mask<Polyhedron_3> mask, int step = 1); 

/*!

applies the PTQ refinement on the control mesh `p` `step` times, 
where `p` contains only triangle facets. 
The geometry of the refined mesh is computed by the geometry policy `mask`. 
This function overwrites the control mesh `p` with the refined mesh. 
The result of a non-triangle mesh `p` is undefined. 
*/ 

template <class Polyhedron_3, template <typename> class Mask> 
void PTQ(Polyhedron_3& p, Mask<Polyhedron_3> mask, int step = 1); 

/*!

applies the DQQ refinement on the control mesh `p` `step` times. 
The geometry of the refined mesh is computed by the geometry policy `mask`. 
This function overwrites the control mesh `p` with the refined mesh. 
*/ 

template <class Polyhedron_3, template <typename> class Mask> 
void DQQ(Polyhedron_3& p, Mask<Polyhedron_3> mask, int step = 1); 

/*!

applies the \f$ \sqrt{3}\f$ triangulation on the control mesh `p` 
`step` times, where `p` contains only triangle facets. 
The geometry of the refined mesh is computed by the geometry policy `mask`. 
This function overwrites the control mesh `p` with the refined mesh. 
The result of a non-triangle mesh `p` is undefined. 
*/ 

template <class Polyhedron_3, template <typename> class Mask> 
void Sqrt3(Polyhedron_3& p, Mask<Polyhedron_3> mask, int step = 1); 

/*!

applies Catmull-Clark subdivision `step` times on the control mesh `p`. 
This function overwrites the control mesh `p` with the subdivided mesh. 
*/ 

template <class Polyhedron_3> 
void CatmullClark_subdivision(Polyhedron_3& p, int step = 1); 

/*!

applies Loop subdivision `step` times on the control mesh `p`. 
This function overwrites the control mesh `p` with the subdivided mesh. 
*/ 

template <class Polyhedron_3> 
void Loop_subdivision(Polyhedron_3& p, int step = 1); 

/*!

applies Doo-Sabin subdivision `step` times on the control mesh `p`. 
This function overwrites the control mesh `p` with the subdivided mesh. 
*/ 

template <class Polyhedron_3> 
void DooSabin_subdivision(Polyhedron_3& p, int step = 1); 

/*!

applies \f$ \sqrt{3}\f$ subdivision `step` times on the control mesh `p`. 
This function overwrites the control mesh `p` with the subdivided mesh. 
*/ 

template <class Polyhedron_3> 
void Sqrt3_subdivision(Polyhedron_3& p, int step = 1); 


/// @}

} /* end namespace Subdivision_method_3 */
} /* end namespace CGAL */

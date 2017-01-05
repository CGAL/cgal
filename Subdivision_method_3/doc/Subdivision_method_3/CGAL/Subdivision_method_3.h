
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


\sa `CGAL::CatmullClark_mask_3<PolygonMesh>`
\sa `CGAL::Loop_mask_3<PolygonMesh`
\sa `CGAL::Sqrt3_mask_3<PolygonMesh>`
*/
/// @{

/*!
 *
 * applies the PQQ refinement several times on the control mesh `pmesh`. 
 * The geometry of the refined mesh is computed by the geometry policy `mask`. 
 * This function overwrites the control mesh `pmesh` with the refined mesh.
 *
 * @tparam PolygonMesh a model of `HalfedgeListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 **/ 

template <class PolygonMesh, class NamedParameters> 
void PQQ(PolygonMesh& pmesh, const NamedParameters& np); 

/*!
 *
 * applies the PTQ refinement several times on the control mesh `pmesh`. 
 * The geometry of the refined mesh is computed by the geometry policy `mask`. 
 * This function overwrites the control mesh `pmesh` with the refined mesh.
 *
 * @tparam PolygonMesh a model of `HalfedgeListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 **/ 

template <class PolygonMesh, class NamedParameters> 
void PTQ(PolygonMesh& pmesh, const NamedParameters& np); 

/*!
 *
 * applies the DQQ refinement several times on the control mesh `pmesh`. 
 * The geometry of the refined mesh is computed by the geometry policy `mask`. 
 * This function overwrites the control mesh `pmesh` with the refined mesh.
 *
 * @tparam PolygonMesh a model of `HalfedgeListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \precondition `pmesh` must be a triangle mesh.
 **/ 

template <class PolygonMesh, class NamedParameters> 
void DQQ(PolygonMesh& pmesh, const NamedParameters& np); 

/*!
 *
 * applies the \f$ \sqrt{3}\f$ refinement several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy `mask`. 
 * This function overwrites the control mesh `pmesh` with the refined mesh. 
 *
 * @tparam PolygonMesh a model of `HalfedgeListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \precondition `pmesh` must be a triangle mesh.
 **/ 

template <class PolygonMesh, class NamedParameters> 
void Sqrt(PolygonMesh& pmesh, const NamedParameters& np); 




/*!
 *
 * applies Catmull-Clark subdivision several times on the control mesh `pmesh`. 
 * This function overwrites the control mesh `pmesh` with the subdivided mesh. 
 
 * @tparam PolygonMesh a model of `HalfedgeListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \precondition `pmesh` must be a triangle mesh.
 **/ 

  template <class PolygonMesh, class NamedParameters> 
void CatmullClark_subdivision(PolygonMesh& pmesh, const NamedParameters& np); 

/*!
 *
 * applies Loop subdivision several times on the control mesh `pmesh`. 
 * This function overwrites the control mesh `pmesh` with the subdivided mesh. 
 
 * @tparam PolygonMesh a model of `HalfedgeListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 **/ 

  template <class PolygonMesh, class NamedParameters> 
void Loop_subdivision(PolygonMesh& pmesh, const NamedParameters& np); 


/*!
 *
 * applies DooSabin subdivision several times on the control mesh `pmesh`. 
 * This function overwrites the control mesh `pmesh` with the subdivided mesh. 
 
 * @tparam PolygonMesh a model of `HalfedgeListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 **/ 

  template <class PolygonMesh, class NamedParameters> 
void DooSabin_subdivision(PolygonMesh& pmesh, const NamedParameters& np); 

 
/*!
 *
 * applies \f$ \sqrt{3}\f$-subdivision several times on the control mesh `pmesh`. 
 * This function overwrites the control mesh `pmesh` with the subdivided mesh. 
 *
 * @tparam PolygonMesh a model of `HalfedgeListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \precondition `pmesh` must be a triangle mesh.
 **/ 

template <class PolygonMesh, class NamedParameters> 
void Sqrt3_subdivision(PolygonMesh& pmesh, const NamedParameters& np); 


/// @}

} /* end namespace Subdivision_method_3 */
} /* end namespace CGAL */
